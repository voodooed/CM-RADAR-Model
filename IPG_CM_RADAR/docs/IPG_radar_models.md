# The two CarMaker radar models

**Start here.** It builds up what a radar actually does, then shows why IPG ships two
completely different radar models and what each one is really computing.

Every worked number below comes from CarMaker's own shipped example sensors and
has been checked against the reference implementation in this package.

Reading time: ~25 minutes. Where to go next is at the end.

---

## Part I — What a radar does

### 1. The basic trick

A radar transmits a radio wave, the wave bounces off something, and a little of
it comes back. Everything a radar knows comes from comparing what came back
with what went out. Three comparisons matter:

| Compare | Gives you | Why |
|---|---|---|
| **when** it came back | **range** | the wave travels at the speed of light |
| **how its frequency shifted** | **radial velocity** | Doppler effect |
| **which antenna element saw it first** | **angle** | the wave arrives tilted across the antenna |

Automotive radar uses 77 GHz, so the wavelength is

```
lambda = c / f = 299 792 458 / 77e9 = 3.89 mm
```

Remember that number — 3.89 mm is roughly a grain of rice, and it explains most
of radar's odd behaviour.

### 2. Range

The wave goes out and comes back, so it covers twice the distance:

```
R = c * tau / 2
```

A target at 100 m returns its echo **667 nanoseconds** later. You cannot measure
that with a stopwatch, so real radars sweep the transmit frequency over time (a
"chirp", FMCW) and measure the *frequency difference* between what is being
transmitted now and the delayed echo. A frequency difference is easy to measure.

The consequence you need to remember: **range resolution is set by how much
bandwidth you sweep.**

```
dR = c / (2 B)
```

* CarMaker's example object-list radar declares `ResolutionDistance = 1.8 m` —
  that corresponds to about 83 MHz of sweep.
* CarMaker's example RSI radar has `NoiseBandwidth = 200 MHz` → 0.75 m.

Two targets closer together than `dR` merge into one. This is not a software
limitation; it is physics, and **both** CarMaker models reproduce it (in very
different ways, as we'll see).

### 3. Velocity

If the target moves, the echo comes back at a slightly different frequency:

```
f_doppler = 2 * v_radial / lambda
```

At 77 GHz, **1 m/s of closing speed shifts the echo by 514 Hz**; 30 m/s shifts
it by 15.4 kHz. That is a big, easily measured shift — which is why radar is so
good at velocity and why it separates a moving car from a stationary guardrail
so easily.

Note the word **radial**. Radar only measures the component of velocity along
the line of sight. A car crossing exactly perpendicular in front of you has
*zero* Doppler and looks, to the Doppler processor, exactly like a parked car.

### 4. Angle

A single antenna cannot tell direction. Two tricks are used:

* **A narrow beam.** Make the antenna big compared to the wavelength and it only
  illuminates a narrow cone. Beamwidth ≈ `lambda / aperture_size`. For a 20°
  beam at 77 GHz you need an aperture of only **9.9 mm** — this is why
  automotive radar antennas are tiny.
* **An array.** Put several receivers side by side. A wave arriving at an angle
  reaches them at slightly different times, so their signals differ in *phase*.
  From those phase differences you can compute the angle. With `N` receivers
  spaced `lambda/2` apart the angular resolution is roughly `2/N` radians:
  16 receivers → 7.2°, 60 receivers → 1.9°.

**This is the crucial asymmetry of radar**: range and velocity are measured
superbly (centimetres, centimetres per second), angle is measured poorly
(degrees). Anything that surprises you about radar usually traces back to this.

### 5. How much power comes back — the radar equation

This is the single most important formula in the whole subject. Build it in four
steps.

**Step 1 — the wave spreads out.** Power `P` radiated in all directions, at
distance `R`, is spread over a sphere of area `4 pi R^2`. If the antenna focuses
it into a beam, it is `G` times more concentrated in that direction:

```
power density at the target  =  P * G / (4 pi R^2)
```

**Step 2 — the target intercepts some and re-radiates it.** How much is
described by one number, the **radar cross section** `RCS` (units: m²). By
definition, the target behaves as if it captured `RCS` square metres of the wave
and re-radiated it uniformly in all directions:

```
power density back at the radar  =  [P * G / (4 pi R^2)] * RCS / (4 pi R^2)
```

**Step 3 — the receiving antenna catches some.** Its effective catching area is
`G * lambda^2 / (4 pi)`.

**Step 4 — multiply it all together:**

```
             P * G^2 * lambda^2 * RCS
    S  =  ----------------------------
              (4 pi)^3 * R^4
```

This is **exactly CarMaker's Equation 505** (it adds two loss factors: `L_A` for
system losses, `L_atm` for rain/fog/atmosphere).

The thing to internalise is the **`R^4`**. The wave spreads out on the way there
*and* on the way back. Consequences:

* Doubling the range costs **12 dB** of signal.
* To double your detection range you need **16 times** the transmit power.
* Conversely, 6 dB of extra margin only buys you 41 % more range.

That is why radar detection ranges are so stubborn, and why a small change in
`RCS` matters far less than intuition suggests.

### 6. Noise, and why detection is a *decision*

The receiver is warm, and warm electronics generate noise. The floor is

```
N = k_B * T * B * F
```

`k_B` = Boltzmann's constant, `T` = temperature, `B` = bandwidth, `F` = the
receiver's "noise figure" (how much worse than ideal it is). This is
**CarMaker's Equation 506**.

For CarMaker's example object-list radar (`T = 293 K`, `B = 25 kHz`,
`F = 4.8 dB`), the noise floor is **−155.15 dBW**.

Now put the two together. That radar (14 dBm transmit, 20.4 dB antenna gain)
looking at a 1 m² target at 100 m receives:

```
signal  = -136.44 dBW
noise   = -155.15 dBW
SNR     =   18.71 dB
```

So the echo is about 74 times stronger than the noise. Good.

But the noise is *random*. Sometimes it is unusually large and looks like a
target; sometimes it partially cancels a real echo and the target disappears.
**Detection is therefore not a yes/no fact about the world — it is a statistical
decision.** You set a threshold, and you trade two errors against each other:

* set it low → you detect weak targets, but noise spikes cause **false alarms**
* set it high → few false alarms, but you **miss** weak targets

The two knobs are the **probability of false alarm** `P_FA` and the
**probability of detection** `P_D`. Standard detection theory relates them to
the required SNR; CarMaker's **Equation 504** is exactly that relation:

```
SNR_min = 2 * ( erfc^-1(2 P_FA) - erfc^-1(2 P_Dmin) )^2
```

With CarMaker's defaults (`P_FA = 1e-6`, `P_Dmin = 0.5`) this gives
**`SNR_min` = 13.54 dB**. Our 1 m² target at 100 m had 18.71 dB, so it is
detected — with `P_D = 0.99995`. Push it out to about 134 m and its SNR falls to
13.54 dB, where it is detected only half the time. That is the detection range.

### 7. RCS is not a property of an object

`RCS` deserves a warning. It is *not* "how big the object is". It is how
effectively the object redirects energy **back toward the radar**, and it depends
violently on geometry:

* A flat metal plate seen head-on is a mirror pointing straight at you —
  enormous RCS. Tilt it 10° and the energy goes elsewhere — RCS collapses.
* A corner (two or three perpendicular surfaces) retro-reflects no matter the
  angle — that is why wheel arches, licence plates and truck underbodies are
  radar-bright.
* A pedestrian is small, soft and irregular — low RCS, and it *flickers* as they
  walk.

Two consequences that both CarMaker models take seriously:

1. **RCS depends strongly on aspect angle.** A car viewed from the side has a
   far larger RCS than one viewed from the rear.
2. **RCS fluctuates.** A real target is many small scatterers whose echoes add
   with random relative phases; move a few millimetres (a fraction of 3.89 mm!)
   and they re-add differently. The standard statistical description is the
   **Swerling models**; CarMaker uses **Swerling 1**, which means the RCS is
   drawn from an exponential distribution around its average
   (**Equation 509**).

This is why a real radar's detections blink in and out even for a perfectly
steady target.

---

## Part II — Two ways to simulate all this

Everything in Part I is physics. A simulator now has to choose **at what level to
reproduce it**, and this is where the two CarMaker models part ways.

```
   REALITY                          SIMULATION OPTION A          SIMULATION OPTION B
   -------                          -------------------          -------------------
   wave leaves antenna              (skipped)                    launch rays
   bounces off the world            (skipped)                    trace reflections
   echo returns                     (skipped)                    sum complex fields
   receiver samples it              (skipped)                    build data cube + noise
   signal processing (FFTs)         (skipped)                    run the FFTs
   threshold / CFAR                 (skipped)                    run the CFAR
   detections clustered             (skipped)                    -> detections
   tracker -> object list           "would this OBJECT be        (a tracker would sit
                                     detected? yes/no"            on top, not modelled)
                                            |                            |
                                            v                            v
                                    RADAR SENSOR (Model 1)       RADAR RSI (Model 2)
                                    object list out              point cloud out
```

**Option A** asks one question per object: *given the radar equation and the
noise floor, would this object's echo clear the threshold?* If yes, report the
object with realistic errors. It never simulates a waveform. It is cheap enough
to run in real time on a CPU.

**Option B** simulates the electromagnetic wave in the actual 3D scene, then
runs a real signal-processing chain over the result. Detections *emerge*; nobody
ever asks "is this a car?". It needs a GPU.

Neither is "better". They answer different questions — see §12.

---

## Part III — Model 1: the Radar Sensor (object list)

CarMaker calls it `Type = "Radar"`. It produces the kind of object list an ADAS
function (ACC, AEB) actually consumes.

### 8. The object model

Every object in the world is reduced to two things:

* an **oriented bounding box** (length × width × height), and
* an **RCS table**: a curve of RCS versus the angle you are looking at it from,
  in the object's *own* frame. CarMaker ships measured tables for
  `RCS_Car`, `RCS_Truck`, `RCS_Pedestrian`, `RCS_Bicycle`, `RCS_GuardRailPost`.

That's it. No mesh, no materials. An object with no RCS table is **invisible** to
this model — which is worth knowing when a custom traffic object mysteriously
never gets detected.

The **centre of the bounding box** is the reference point for range, angle and
velocity. (Remember this — Model 2 behaves differently, see §11.)

### 9. What it computes, per object, per cycle

```
 1. Is the object inside the field of view and the range gate?   -> candidate
 2. Look up the antenna gain G in the direction of the object
 3. Look up RCS at the aspect angle, then
      - randomise it (Swerling 1)              <- the target flickers
      - reduce it if something blocks the view <- occlusion
 4. Apply rain / fog / atmospheric loss
 5. Radar equation (Eq. 505)                   -> received signal S
 6. Noise floor (Eq. 506) + road clutter       -> noise N
 7. SNR = S/N. Detected if SNR > SNR_min (Eq. 504).
 8. Merge objects the radar cannot separate    <- resolution
 9. Add measurement noise to range/angle/speed <- accuracy
10. Bookkeeping: confidence, new/measured/lost, list size
11. Optionally invent false positives
12. Delay the whole result by the sensor latency
```

Three of these steps are what make the model *realistic* rather than ideal, and
they are worth dwelling on.

**Step 7 is why false negatives happen.** Because RCS is randomised in step 3, an
object sitting near the threshold is detected in some cycles and not in others.
A pedestrian at 40 m in CarMaker's example configuration has an SNR of about
14 dB against a 13.54 dB threshold — so they flicker, exactly as a real radar
would show them. This is not a bug in the model; it is the point of the model.

**Step 8 is the resolution cell.** Two objects are reported as one if they are
too close in **range and velocity and angle simultaneously**:

```
separable if  |dR| > delta*R_range  OR  |dv| > delta*R_vel  OR  |dangle| > delta*R_angle
```

`delta` is the `Separability` parameter (1.5 in the example). With
`ResolutionDistance = 1.8 m`, two cars less than 2.7 m apart in range — and
also similar in speed and bearing — merge into one object with one bounding box.
This is a classic ADAS failure mode and the model reproduces it.

**Step 10 is the tracker.** The model keeps a confidence counter per object
(`ProbExist`, 0–7 ≈ 0 %…100 %). It goes up while the object is detected and down
when it is not, and the object stays on the list until the counter hits zero.
So an object that flickers out for one cycle does **not** vanish from the output —
which is exactly how a real radar's tracker behaves, and something a naive
"is it visible?" sensor model gets wrong.

### 10. What it deliberately does *not* model

No multipath, no ghost targets from physics, no micro-Doppler, no imaging of
individual scattering centres on a vehicle. False positives exist, but as a
*statistical sprinkle* (a Poisson number of clutter objects, plus mirror objects
near guardrails) rather than as a consequence of simulated wave propagation.

If you need those effects, you need Model 2.

---

## Part IV — Model 2: the Radar RSI (raw signal)

CarMaker calls it `Type = "RadarRSI"`. RSI = Raw Signal Interface.

### 11. What it computes

**Stage 1 — shoot rays.** Tens of thousands of rays are launched across the
field of view (87 230 in CarMaker's shipped example). Each ray carries a
**complex amplitude** and a **polarisation** — complex, because interference is
the whole point.

**Stage 2 — bounce them around the real geometry.** Each ray hits actual triangle
meshes with actual materials (`MaterialLib` gives each material a permittivity
and a roughness). At each hit:

* the reflection strength comes from the **Fresnel equations** and the material's
  permittivity — metal reflects almost everything, asphalt much less;
* the surface normal is randomly perturbed by the material's roughness, which is
  how rough surfaces scatter instead of mirroring;
* the ray continues, and the interaction also sends energy **back toward the
  receiver** — both directly (if there's line of sight) and back along the path
  it came from.

That last point is where **ghost targets** come from. A wave that goes
sensor → guardrail → car → guardrail → sensor arrives back looking as if it came
from the guardrail direction. The radar dutifully reports a car that isn't there.
Model 2 produces those for free, because it is simulating the actual path.

**Stage 3 — accumulate the echoes.** Every interaction point contributes its
complex amplitude at its own range and Doppler. Because the amplitudes are
complex, echoes that arrive in phase **add** and echoes out of phase **cancel**.
This is real: bouncing off the road surface creates a second path a few
millimetres longer than the direct one, and since 3.89 mm is the wavelength,
the two paths swing in and out of phase as the target's range changes. Our
validation measures exactly this — the signal strength oscillating by ±2…4 dB as
a target approaches. Real 77 GHz radars show this and it is a genuine source of
detection dropouts.

**Stage 4 — the device model.** Now the simulation behaves like the actual
electronics:

```
  data cube:  range x Doppler x receiver
       |
       + add thermal noise in every cell (Eq. 525)
       |
  2-D FFT view -> a "Range-Doppler map": a picture where each pixel is
       |          one (range, speed) combination and its brightness is power
       |
  OS-CFAR filter: for each pixel, look at its neighbours, and call it a
       |          detection only if it is much brighter than they are
       |          <- this is how real radars threshold, and it is why a
       |             bright target doesn't create a halo of false detections
       |
  peak finder:   keep only local maxima
       |
  angular FFT across the receivers -> which direction is it in?
       |
  CA-CFAR + peak finder again
       |
  parabolic interpolation -> sub-bin accuracy (Eqs. 527-529)
       v
  a point cloud: (x, y, z) or (range, azimuth, elevation), power, velocity
```

The CFAR ("constant false alarm rate") stage is the Model 2 equivalent of
Model 1's `SNR > SNR_min` — same idea (threshold against noise), but applied per
pixel of an image rather than per object, with the noise level *estimated
locally* from the neighbourhood instead of assumed.

### 12. The one difference that will surprise you

Run the demo in this package and you get, for a 4.5 m car whose **centre** is at
40.6 m:

```
Model 2 reports:  range = 38.45 m
```

That is not an error. 38.35 m is where the car's **front bumper** is. A
ray-traced radar detects surfaces, and the first surface is the bumper. Model 1,
by contrast, always reports the bounding-box centre at 40.6 m.

Generalise this: **Model 1 reports objects; Model 2 reports surfaces.** A single
truck in Model 2 produces a spray of detections along its side; in Model 1 it
produces one entry. If you are developing a clustering or tracking algorithm,
that difference is the entire reason Model 2 exists.

---

## Part V — Choosing, and where to read next

### 13. Which model for which job

| You are working on | Use | Because |
|---|---|---|
| ACC / AEB / lane-change logic | **Model 1** | you need a plausible object list, fast, in real time |
| Sensor-fusion logic consuming object lists | **Model 1** | that's the interface it will see in the car |
| Large scenario campaigns, thousands of runs | **Model 1** | CPU, real-time, cheap |
| Clustering / tracking algorithms | **Model 2** | you need raw detections, not pre-formed objects |
| Radar signal processing, CFAR tuning, MIMO/beamforming | **Model 2** | you need the data cube and the virtual receivers |
| Ghost targets, tunnels, guardrails, NLOS detection | **Model 2** | Model 1 doesn't simulate wave paths |
| Sensor placement / mounting-angle studies | **Model 2** | multipath and occlusion depend on real geometry |
| Understanding a specific real-radar dropout | **Model 2** | usually multipath or a resolution-cell effect |

A common pattern in industry: develop the function against Model 1 because you
can run 10 000 scenarios overnight, then re-run the interesting failures through
Model 2 to see whether the physics would have made them worse.

### 14. The five things worth remembering

1. **`R^4`.** Everything about detection range is dominated by it.
2. **Detection is statistical**, not geometric. Both models threshold an SNR.
3. **Angle is the weak measurement**; range and velocity are excellent.
4. **RCS fluctuates**, so detections flicker — by design, in both models.
5. **Model 1 reports objects at their centre; Model 2 reports surfaces.**

### 15. Where to go next

```
you are here ->  docs/IPG_radar_models.md
                          |
        +-----------------+------------------+
        v                                    v
 radar_object_list/explanation.md     radar_rsi/explanation.md
   what Model 1 does, with              what Model 2 does, with
   the actual IPG equations             the actual IPG equations
        |                                    |
        v                                    v
 radar_object_list/                   radar_rsi/
   implementation_logic.md              implementation_logic.md
   how to build it                      how to build it
        |                                    |
        +-----------------+------------------+
                          v
             validation/validation_report.md
             what we could prove, and what we had to assume
```

Then run the two demos and read their output next to this document:

```bash
# from the package root; CARMAKER_DIR is optional (see README -> Running it)
export CARMAKER_DIR=<path-to>/carmaker/<arch>-<version>
python3 radar_object_list/example.py --ipg "$CARMAKER_DIR"
python3 radar_rsi/example.py --ipg "$CARMAKER_DIR" --rays 4000
```

**IPG Reference Manual** Use as a lookup. Its radar chapters
(*Sensors → Radar Sensor* and *Sensors → Radar RSI*) are short and assume the
background in Part I. The two `explanation.md` files quote every relevant
equation and cite the exact section, so go to the manual only when you want to
confirm a specific sentence — `docs/source_inventory.md` lists the file path of
every section used.

**One caution for all of this package.** CarMaker's radar models ship as
compiled code, and IPG documents *what* they compute far more completely than
*how* in a few places (road clutter, atmospheric damping coefficients, the exact
scattering kernel). Everything reconstructed rather than documented is marked
`[INFERRED]` in the code and listed explicitly in the two `explanation.md` files.
Trust the documented equations; treat the inferred parts as a reasonable
reconstruction, not as IPG's answer.

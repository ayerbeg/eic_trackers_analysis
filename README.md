# lambda_analysis

ROOT macro for $\Lambda \to p + \pi^-$ decay vertex reconstruction in the EIC far-forward region.

## Physics goal

Establish the minimum spatial resolution requirements for a charged-particle tracking detector placed in the field-free region between the last beamline magnet (~20.8 m) and the Zero-Degree Calorimeter (ZDC) at the EIC. The tracker is intended to reconstruct the decay vertex of $\Lambda \to p + \pi^-$ from the two daughter tracks, enabling the statement:

> *"Two tracks originate from a common point at z = Z mm, consistent with a Λ decay in the field-free region."*

No momentum measurement is available in the field-free region — invariant mass reconstruction is outside scope.

## Analysis summary

- **Beam energy:** 18×275 GeV (electron × hadron)
- **Simulation:** dd4hep/GEANT4 via EIC simulation framework
- **Data format:** EDM4hep/PODIO ROOT files
- **Signal:** Λ → p + π⁻ (BR ≈ 64%), primary Lambda from hard scatter
- **Acceptance:** Λ decaying in the field-free region, 20.8 m < z < 35 m
- **Tracker geometry:** 4 planes over 2 m, 500 mm upstream of ZDC face
- **Pitch scenarios:** 50 µm, 100 µm, 250 µm, 500 µm, 1 mm

## Dependencies

| Library | Purpose |
|---|---|
| ROOT | Histograms, canvases, file I/O |
| PODIO | Reading EDM4hep event frames |
| EDM4hep | MC particle data model |
| fmt | Formatted stdout output |

All dependencies are available in the standard EIC software environment.

## How to run

```bash
# Full run (all events, all files)
root -x -l -b -q lambda_analysis.cxx

# Interactive (with ROOT GUI)
root -x -l lambda_analysis.cxx
root [0] lambda_analysis()

# Quick test (100 events per file)
root [0] lambda_analysis(100)
```

## Output files

| File | Contents |
|---|---|
| `lambda_analysis_18x275.root` | All histograms, reopenable for replotting |
| `lambda_analysis_18x275.pdf` | Lambda kinematics and daughter distributions |
| `tracker_study_18x275.pdf` | Tracker resolution study (12 pages) |

## Tracker study PDF — page guide

| Page | Content |
|---|---|
| 1 | ΔR at each tracker plane (truth) |
| 2 | ΔR zoomed 0–3 mm — sub-mm region where pitch matters |
| 3–4 | σ_z and σ_α distributions (analytical) |
| 5–6 | σ_z vs α and σ_z vs D_eff (2D, analytical) |
| 7 | **Money plot:** mean σ_z and analytical efficiency vs pitch |
| 8 | **Money plot:** z_reco distributions + smearing efficiency vs pitch |
| 9–10 | Δz and Δα residuals (smearing) |
| 11 | Δz vs α (2D, smearing) |
| 12 | Validation: RMS(Δz) vs analytical σ_z per α bin |

## Configuration

All tunable parameters are in two namespaces at the top of the file. **No other part of the code needs to be changed for typical modifications.**

### `namespace Cuts` — physics selection

```cpp
constexpr int    LAMBDA_PRIMARY_IDX = 4;      // generator-dependent index
constexpr double DECAY_Z_MIN        = 20800;  // mm — last magnet exit
constexpr double DECAY_Z_MAX        = 35000;  // mm — ZDC face
constexpr double CROSSING_ANGLE     = 0.025;  // rad — 25 mrad
constexpr int    N_FILES            = 20;
```

### `namespace Tracker` — detector geometry

```cpp
constexpr int    N_PLANES       = 4;
constexpr double LENGTH         = 2000;   // mm
constexpr double STANDOFF       = 500;    // mm from ZDC
constexpr double SIGMA_X[]      = {0.05, 0.10, 0.25, 0.50, 1.00};  // mm
constexpr double EFF_THRESH[]   = {500, 1000, 2000};  // mm
constexpr int    PRIMARY_THRESH = 1;      // index — 1000 mm is primary
```

## Key results (18×275 GeV, 20000 events)

**Lambda selection:**
- Total events: 20,000
- Decayed in field-free region (all channels): 2,211
  - Charged p+π⁻ (used in analysis): 1,324
  - Neutral n+π⁰ (not analysed): 721
  - Material interactions: 166

**Track separation:** ΔR at tracker planes peaks at ~20–25 mm, far larger than any pitch scenario. Track separation is never the limiting factor — all scenarios resolve both tracks at plane 1.

**Efficiency comparison at 1 m threshold:**

| Pitch | Analytical | Smearing |
|---|---|---|
| 50 µm | 99.4% | 94.6% |
| 100 µm | 97.3% | 87.3% |
| 250 µm | 87.0% | 71.8% |
| 500 µm | 52.2% | 52.7% |
| 1 mm | 24.1% | 34.5% |

The analytical formula is optimistic for fine pitch (50–250 µm) and pessimistic for coarse pitch (500 µm–1 mm). **The smearing efficiency is the more honest design metric.** To achieve 80% smearing efficiency at the 1 m threshold, a pitch of ~100–150 µm is required.

## How to add a new histogram

1. Book it in `book_histograms()`:
   ```cpp
   H["MyHist"] = new TH1F("MyHist", "Title;X;Y", 100, 0, 100);
   ```
2. Fill it in `process_event()`:
   ```cpp
   h1("MyHist")->Fill(my_value);
   ```
3. Draw it in `plot()` or `plot_tracker_study()`.
4. It is automatically saved to the ROOT file — no extra step needed.

## How to add a new pitch scenario

1. Add the value to `Tracker::SIGMA_X[]`
2. Increment `Tracker::N_SC`
3. Add a label to `Tracker::LABELS[]`

All loops update automatically.

## Documentation

| Document | Contents |
|---|---|
| `lambda_analysis_code_reference.pdf` | Full code reference manual |
| `tracker_resolution_note_v2.pdf` | Physics technical note with derivations |
| `derivations.pdf` | Mathematical derivations (σ_θ, σ_z, α) |
| `references.bib` | BibTeX references |

## Author

C. Ayerbe Gayoso — EIC Meson Structure Group  
Analysis for the EIC far-forward detector instrumentation programme.

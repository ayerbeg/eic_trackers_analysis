// ============================================================================
// lambda_analysis.cxx
//
// Analysis of Lambda (PDG 3122) -> p + pi- decays in the EIC far-forward
// region, between the last beamline magnet and the Zero Degree Calorimeter
// (ZDC). Reads EDM4hep MC truth from PODIO ROOT files produced by the EIC
// simulation framework (dd4hep/GEANT4).
//
// Physics goal: estimate the spatial resolution needed in a tracker placed
// just in front of the ZDC to reconstruct the p and pi- from Lambda decay.
//
// Entry point for ROOT:
//   root -x -l lambda_analysis.cxx
//   root -x -l -b -q lambda_analysis.cxx          // batch, no GUI
//   root [0] lambda_analysis(100)                  // 100 events per file
// 
// VERSION HISTORY:
//
// - Initial version: Dmitry Romanov from the EIC Meson Structure Working Group
// https://github.com/JeffersonLab/meson-structure/blob/main/tutorials/cpp_01_read_edm4eic.cxx
// 
// - 2026-03 Modified by C. Ayerbe to analyze Lambda decays in the far-forward region and 
// add new histograms for the daughter kinematics at the ZDC.
// - 2026-05 Further edits by C. Ayerbe using CLAUDE to clean up the code structure, add
// comments, and implement the new DaughterKin struct and Counters struct for better organization.
// - 2026-06 Further edits by C. Ayerbe to implement the tracker resolution study, including new histograms 
// and the readplot_tracker() function for the ROOT macro entry point.
// - 2026-06 C. Ayerbe + Claude: added Tracker namespace (4 planes, 2 m, 500 mm standoff,
// 5 sigma_x scenarios); analytical sigma_z and sigma_alpha per event; Gaussian smearing
// with PCA vertex reconstruction; efficiency calculation (fraction of events with
// sigma_z < 500/1000/2000 mm); plot_tracker_study() and save_tracker_study() producing
// a separate 8-page PDF.
// - 2026-06 C. Ayerbe + Claude: added TProfile ResZ_rms_vsAlpha per scenario
// (option "S" stores RMS of Delta_z in alpha slices); added validation canvas
// (page 9) overlaying analytical sigma_z profile vs smearing RMS per alpha bin
// to validate the Gaussian approximation and confirm the formula.
// - 2026-06 C. Ayerbe + Claude: added LambdaEnd_analysis histogram filled only after all
// analysis cuts (primary, field-free, charged channel); expanded Counters struct with
// per-rejection-reason breakdown (lambda_escaped, lambda_before_magnet, lambda_beyond_zdc,
// lambda_field_free, lambda_charged, lambda_neutral, lambda_other_ndau); neutral channel
// daughters now filled in a dedicated block before the charged-channel cut; diagnostic
// canvas expanded to 2x3 with overlay of field-free vs analysis sample (pad 5) and
// clean analysis-sample distribution (pad 6).
// ============================================================================

#ifdef __CLING__
R__LOAD_LIBRARY(podioDict)
R__LOAD_LIBRARY(podioRootIO)
R__LOAD_LIBRARY(libedm4hepDict)
R__LOAD_LIBRARY(libedm4eicDict)
#endif

#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include <edm4hep/MCParticleCollection.h>
#include <fmt/core.h>
#include <fmt/format.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TRandom.h>
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace std;

// ============================================================================
// Analysis cuts
// All geometric thresholds and PDG codes live here.
// To change a cut, edit this block only.
// ============================================================================
namespace Cuts {

    // PDG code of the Lambda baryon
    constexpr int    LAMBDA_PDG         = 3122;

    // Object index of the primary Lambda in the MC particle list.
    // In the current generator settings the primary Lambda always has index 4.
    // This is generator-dependent and may need updating if the config changes.
    constexpr int    LAMBDA_PRIMARY_IDX = 4;

    // Number of daughters for the charged decay Lambda -> p + pi-
    constexpr int    LAMBDA_N_DAUGHTERS = 2;

    // Longitudinal decay vertex window [mm].
    // DECAY_Z_MIN: last beamline magnet exit (~20.8 m). Only Lambdas decaying
    //   beyond this point travel in a field-free region where daughters move
    //   in straight lines — the assumption used for the ZDC projection.
    // DECAY_Z_MAX: ZDC face (~35 m). Lambdas beyond this are not useful.
    constexpr double DECAY_Z_MIN        = 20800.0;  // mm
    constexpr double DECAY_Z_MAX        = 35000.0;  // mm

    // ZDC face position [mm]. Used to compute the propagation distance from
    // the Lambda decay vertex to the detector plane.
    constexpr double ZDC_Z              = 35000.0;  // mm

    // EIC beam crossing angle [radians].
    // The hadron beam crosses the electron beam axis at 25 mrad in the x-z
    // plane toward -x. Used to compute the beam-centred coordinate frame.
    constexpr double CROSSING_ANGLE     = 25.0e-3;  // rad

    // PDG codes of the Lambda decay daughters
    constexpr int    PDG_PROTON         =  2212;
    constexpr int    PDG_PION_MINUS     =  -211;
    constexpr int    PDG_NEUTRON        =  2112;
    constexpr int    PDG_PION_ZERO      =   111;

    // Number of input ROOT files to process
    constexpr int    N_FILES            = 20;

    // File path pattern for fmt::format. Must be constexpr const char*
    // (not std::string) because fmt::format requires a compile-time format
    // string in fmt >= 10.
    constexpr const char* FILE_PATTERN =
        "./data-2026-04/dd4hep/18x275/k_lambda_18x275_5000evt_{:04d}.edm4hep.root";

    // Output PDF base name (without extension)
    constexpr const char* OUTPUT_PDF = "lambda_analysis_18x275";
}

// ============================================================================
// Tracker resolution study — configuration
// ============================================================================
namespace Tracker {
    constexpr int    N_PLANES   = 4;        // number of tracking planes
    constexpr double LENGTH     = 2000.0;   // mm — total length spanned by the planes
    constexpr double STANDOFF   = 500.0;    // mm — gap from ZDC face to the LAST plane
    // Layout along z:
    //   |decay|...|plane1|---LENGTH---|plane4|---STANDOFF---|ZDC|
    // plane1 at z = ZDC_Z - STANDOFF - LENGTH = 32500 mm
    // plane4 at z = ZDC_Z - STANDOFF          = 34500 mm
    // D (extrapolation distance) = dist2ZDC - STANDOFF - LENGTH

    // Spatial resolution scenarios [mm]: 50µm, 100µm, 250µm, 500µm, 1mm
    constexpr double SIGMA_X[]  = {0.05, 0.10, 0.25, 0.50, 1.00};
    constexpr int    N_SC       = 5;

    constexpr const char* LABELS[] = {
        "50 #mum", "100 #mum", "250 #mum", "500 #mum", "1 mm"
    };

    // Efficiency thresholds for sigma_z [mm].
    // 1000 mm is the PRIMARY threshold — sufficient to confirm the decay
    // is in the field-free region with ~1 m precision, which is the
    // physically meaningful requirement for a 14 m decay region.
    // 500 mm is "ambitious" (sub-50cm precision, likely overkill).
    // 2000 mm is "minimum acceptable" (2 m precision, very loose).
    constexpr double EFF_THRESH[]  = {500.0, 1000.0, 2000.0};
    constexpr int    N_THRESH      = 3;
    constexpr int    PRIMARY_THRESH = 1; // index into EFF_THRESH[] = 1000 mm

    constexpr const char* OUTPUT_PDF = "tracker_study_18x275";
}

// Angular resolution [rad] for a straight-line fit through N equally-spaced
// planes over total length L with per-plane spatial resolution sigma_x [mm].
//   sigma_theta = (sigma_x / L) * sqrt(12*(N-1) / (N*(N+1)))
// For N=4: factor = sqrt(12*3/(4*5)) = sqrt(1.8) = 1.34164...
// NOTE: this function assumes all 4 planes resolve the tracks (best case).
// Use sigma_theta_eff() for the physically correct per-event calculation.
inline double sigma_theta_rad(double sigma_x_mm) {
    constexpr double factor = 1.3416407864998738; // sqrt(1.8), exact for N=4
    return (sigma_x_mm / Tracker::LENGTH) * factor;
}

// z position [mm] of tracker plane ip (1-indexed, 1=closest to decay, N=closest to ZDC).
// Planes are equally spaced between plane1 and planeN:
//   plane1: z = ZDC_Z - STANDOFF - LENGTH  (32500 mm)
//   planeN: z = ZDC_Z - STANDOFF           (34500 mm)
inline double plane_z(int ip) {
    double z1 = Cuts::ZDC_Z - Tracker::STANDOFF - Tracker::LENGTH;
    double z4 = Cuts::ZDC_Z - Tracker::STANDOFF;
    return z1 + (ip - 1) * (z4 - z1) / (Tracker::N_PLANES - 1);
}

// Angular resolution [rad] using only planes from first_plane to N_PLANES.
// This is the physically correct formula when the two tracks are not resolved
// at the earlier planes — only the resolved planes contribute to the slope fit.
//
//   N_eff = N_PLANES - first_plane + 1  (number of resolved planes)
//   L_eff = plane_z(N_PLANES) - plane_z(first_plane)  (lever arm)
//
// When N_eff < 2 (only one plane resolves the tracks), no slope can be fitted
// and a large sentinel value (1e9) is returned — event is unresolvable.
// This is conservative: merged hits at earlier planes still constrain the
// Lambda direction but not individual track slopes.
inline double sigma_theta_eff(double sigma_x_mm, int first_plane) {
    int    N_eff = Tracker::N_PLANES - first_plane + 1;
    if (N_eff < 2) return 1e9;
    double L_eff = plane_z(Tracker::N_PLANES) - plane_z(first_plane);
    if (L_eff < 1.0) return 1e9;
    double factor = std::sqrt(12.0 * (N_eff - 1)
                              / ((double)N_eff * (N_eff + 1)));
    return (sigma_x_mm / L_eff) * factor;
}

// ============================================================================
// DaughterKin
//
// Holds kinematic quantities for one decay daughter projected onto the ZDC
// plane via straight-line propagation (valid in the field-free region).
//
// Two coordinate frames are computed:
//
//  LAB FRAME  (xAtZDC, yAtZDC, radius)
//    Full position in the global detector frame. The hadron beam axis is
//    offset by ~-875 mm in x at the ZDC due to the 25 mrad crossing angle.
//    Use this to determine where to physically place the tracker.
//
//  BEAM-CENTRED FRAME  (xBC, yBC, radiusBC)
//    Position relative to the hadron beam axis at the ZDC face. The beam
//    axis offset (= -zv * tan(25 mrad)) is subtracted event by event, so
//    the distribution is centred near (0,0). Use this to characterise the
//    angular spread seen by the detector independently of its placement.
//
// Note: DeltaR between p and pi- is identical in both frames because the
// common vertex offset cancels in the difference (x_p - x_pi).
// ============================================================================
struct DaughterKin {
    double theta    = 0;  // polar angle w.r.t. z-axis [degrees]
    double phi      = 0;  // azimuthal angle in x-y plane [degrees]

    // Lab frame
    double xAtZDC   = 0;  // x at ZDC plane, global detector frame [mm]
    double yAtZDC   = 0;  // y at ZDC plane, global detector frame [mm]
    double radius   = 0;  // sqrt(x^2+y^2) in lab frame [mm]

    // Beam-centred frame
    double xBC      = 0;  // x relative to hadron beam axis at ZDC [mm]
    double yBC      = 0;  // y relative to hadron beam axis at ZDC [mm]
    double radiusBC = 0;  // transverse distance from beam axis [mm]

    bool isSet() const { return theta != 0; }

    // Compute both frames from momentum (px,py,pz), decay vertex (xv,yv,zv),
    // ZDC face position zZDC, and crossing angle alpha [radians].
    void compute(double px, double py, double pz,
                 double xv, double yv, double zv,
                 double zZDC, double alpha) {
        double dz = zZDC - zv;
        double p  = std::sqrt(px*px + py*py + pz*pz);

        theta  = std::acos(pz / p) * 180.0 / M_PI;
        phi    = std::atan2(py, px) * 180.0 / M_PI;

        // Lab frame: full straight-line propagation including vertex offset
        xAtZDC = xv + dz * (px / pz);
        yAtZDC = yv + dz * (py / pz);
        radius = std::sqrt(xAtZDC*xAtZDC + yAtZDC*yAtZDC);

        // Beam-centred frame: subtract the hadron beam axis position at ZDC.
        // The beam axis in x at longitudinal position z is: x_beam = -z*tan(alpha).
        // We evaluate it at zZDC (the detector plane), not at zv, because we
        // want positions relative to where the beam axis crosses the detector.
        double xBeamAtZDC = -zZDC * std::tan(alpha);  // beam axis x at ZDC [mm]
        xBC      = xAtZDC - xBeamAtZDC;
        yBC      = yAtZDC;                             // beam axis y = 0 by definition
        radiusBC = std::sqrt(xBC*xBC + yBC*yBC);
    }
};

// ============================================================================
// Counters — tracks particles passing the analysis selections
// ============================================================================
struct Counters {
    int total_events      = 0;

    // ── Primary Lambda breakdown ───────────────────────────────────────────
    // These explain exactly how many primary Lambdas reach each stage and
    // why the remainder are not used in the tracking analysis.
    int lambda_escaped       = 0; // no daughters (Lambda left the simulation volume)
    int lambda_before_magnet = 0; // decays before DECAY_Z_MIN (inside magnet lattice)
    int lambda_beyond_zdc    = 0; // decays past DECAY_Z_MAX (beyond ZDC face)
    int lambda_field_free    = 0; // decays in field-free window (all channels)
    int lambda_charged       = 0; // field-free + exactly 2 daughters = p + pi-  [USED]
    int lambda_neutral       = 0; // field-free + neutral channel (n + pi0)       [not used]
    int lambda_other_ndau    = 0; // field-free + unexpected daughter count (material interactions)

    // ── Daughter counters (for the charged channel used in analysis) ───────
    int proton       = 0;
    int pion_minus   = 0;
    int neutron      = 0;   // neutral channel daughters (informational)
    int pion_zero    = 0;

    void print() const {
        fmt::print("\n========== Lambda Selection Summary ==========\n");
        fmt::print("Total events processed           : {:6d}\n", total_events);
        fmt::print("\n  Primary Lambda disposition:\n");
        fmt::print("    Escaped (no daughters)       : {:6d}\n", lambda_escaped);
        fmt::print("    Decayed before magnet exit   : {:6d}  (z < {:.0f} mm)\n",
                   lambda_before_magnet, Cuts::DECAY_Z_MIN);
        fmt::print("    Decayed beyond ZDC           : {:6d}  (z > {:.0f} mm)\n",
                   lambda_beyond_zdc, Cuts::DECAY_Z_MAX);
        fmt::print("    Decayed in field-free region : {:6d}  (all channels)\n",
                   lambda_field_free);
        fmt::print("      -> charged (p + pi-)       : {:6d}  [used in analysis]\n",
                   lambda_charged);
        fmt::print("      -> neutral (n + pi0)       : {:6d}  [not analysed]\n",
                   lambda_neutral);
        fmt::print("      -> other daughter count    : {:6d}  [material interactions, excluded]\n",
                   lambda_other_ndau);
        fmt::print("\n  Daughter counts (charged channel):\n");
        fmt::print("    Protons found                : {:6d}\n", proton);
        fmt::print("    Pi- found                    : {:6d}\n", pion_minus);
        fmt::print("    Neutrons (neutral channel)   : {:6d}\n", neutron);
        fmt::print("    Pi0 (neutral channel)        : {:6d}\n", pion_zero);
        fmt::print("==============================================\n\n");
    }
} counters;

// ============================================================================
// Histogram registry
//
// All histograms live in one map keyed by a short name string.
// To add a histogram:
//   1. One new TH1F/TH2F line in book_histograms()
//   2. One Fill() call in process_event()
//   3. One Draw() call in plot()
// Nothing else needs to change.
// h1() / h2() are typed accessors that avoid casting at every fill site.
// ============================================================================
map<string, TH1*> H;   // TH1* is the common base class of TH1F and TH2F

// Tracks the file currently being processed — set in process_file() so
// that process_event() can include it in diagnostic printouts.
string current_file = "";

inline TH1F* h1(const string& name) { return static_cast<TH1F*>(H.at(name)); }
inline TH2F* h2(const string& name) { return static_cast<TH2F*>(H.at(name)); }

void book_histograms() {

    // ── Lambda decay endpoint along z ───────────────────────────────────────
    // Where the primary Lambda stops existing (i.e. where it decays).
    // Sub-ranges: inside the magnet lattice (0-6 m), between magnets
    // (6-20.8 m), and in the field-free region (>20.8 m).
    H["LambdaEnd"]   = new TH1F("LambdaEnd",
        "Lambda decay endpoint (all);z_{end} (mm);Count",      1000, 0, 35000);
    H["LambdaEnd_1"] = new TH1F("LambdaEnd_1",
        "Lambda decay endpoint (0-6 m);z_{end} (mm);Count",    1000, 0, 6000);
    H["LambdaEnd_2"] = new TH1F("LambdaEnd_2",
        "Lambda decay endpoint (6-20.8 m);z_{end} (mm);Count", 1000, 6000, 20800);
    H["LambdaEnd_3"] = new TH1F("LambdaEnd_3",
        "Lambda decay endpoint (primary #Lambda, all channels, field-free)"
        ";z_{v} (mm);Count",
        200, Cuts::DECAY_Z_MIN, Cuts::DECAY_Z_MAX);

    // ── Lambda decay endpoint: analysis sample only ─────────────────────────
    // Filled ONLY for primary Lambdas that pass ALL analysis cuts:
    //   - primary index = 4
    //   - decay in (DECAY_Z_MIN, DECAY_Z_MAX)
    //   - exactly 2 daughters (charged channel only)
    // The integral of this histogram equals counters.lambda_charged.
    // Compare with LambdaEnd_3 to see what fraction of field-free Lambdas
    // are actually used (charged channel fraction and z-window efficiency).
    H["LambdaEnd_analysis"] = new TH1F("LambdaEnd_analysis",
        "Lambda decay endpoint [analysis sample: primary, charged, field-free]"
        ";z_{v} (mm);Count",
        200, Cuts::DECAY_Z_MIN, Cuts::DECAY_Z_MAX);

    // ── Lambda production vertex along z ────────────────────────────────────
    // Where the Lambda was created. For primary Lambdas from the hard scatter
    // this should peak near z=0 (the interaction point).
    H["LambdaGenZ"]  = new TH1F("LambdaGenZ",
        "Lambda gen vertex z;z_{vtx} (mm);Count",              1000, 0, 35000);
    H["LambdaGen_1"] = new TH1F("LambdaGen_1",
        "Lambda gen vertex (0-6 m);z_{vtx} (mm);Count",        100, 0, 6000);
    H["LambdaGen_2"] = new TH1F("LambdaGen_2",
        "Lambda gen vertex (6-20.8 m);z_{vtx} (mm);Count",     100, 6000, 20800);
    H["LambdaGen_3"] = new TH1F("LambdaGen_3",
        "Lambda gen vertex (>20.8 m);z_{vtx} (mm);Count",      100, 20800, 35000);

    // ── Lambdas that never decayed inside the simulation volume ─────────────
    // getDaughters().size()==0 means the Lambda escaped before decaying.
    // Their endpoint distribution shows how many are "lost" beyond the ZDC.
    H["LambdaZnoDaughters"] = new TH1F("LambdaZnoDaughters",
        "Lambda endpoint (no daughters);z_{end} (mm);Count",   1000, 0, 35000);

    // ── Lambda polar angle at production ────────────────────────────────────
    // Small angles (~mrad) expected for forward-produced Lambdas at EIC.
    H["LambdaAngleGen"] = new TH1F("LambdaAngleGen",
        "Lambda polar angle;#theta (deg);Count",                100, 0, 10);

    // ── Daughter: proton (PDG 2212) ─────────────────────────────────────────
    H["Vertex_p"]    = new TH1F("Vertex_p",
        "Proton vertex z;z_{vtx} (mm);Count",                  1000, 0, 35000);
    H["Endpoint_p"]  = new TH1F("Endpoint_p",
        "Proton endpoint z;z_{end} (mm);Count",                1000, 0, 35000);
    H["Theta_p"]     = new TH1F("Theta_p",
        "Proton polar angle;#theta (deg);Count",                100, 0, 10);
    // Lab frame: centred near x ~ -875 mm due to 25 mrad crossing angle
    H["XYatZDC_p"]   = new TH2F("XYatZDC_p",
        "Projected XY at ZDC (p) [lab];X (mm);Y (mm)",         100,-1200,-600, 100,-300, 300);
    // Beam-centred: beam axis subtracted, centred near (0,0)
    H["XYatZDC_p_BC"] = new TH2F("XYatZDC_p_BC",
        "Projected XY at ZDC (p) [beam-centred];X_{BC} (mm);Y_{BC} (mm)", 100,-300,300, 100,-300,300);
    H["ZvsTheta_p"]  = new TH2F("ZvsTheta_p",
        "Z from ZDC vs #theta (p);#theta (deg);Z from ZDC (mm)", 100,0,10, 100,0,15000);
    H["ZvsRadius_p"] = new TH2F("ZvsRadius_p",
        "Z from ZDC vs R at ZDC (p);R_{ZDC} (mm);Z from ZDC (mm)", 100,0,1500, 100,0,15000);

    // ── Daughter: pi- (PDG -211) ────────────────────────────────────────────
    H["Vertex_pi"]    = new TH1F("Vertex_pi",
        "#pi^{-} vertex z;z_{vtx} (mm);Count",                 1000, 0, 35000);
    H["Endpoint_pi"]  = new TH1F("Endpoint_pi",
        "#pi^{-} endpoint z;z_{end} (mm);Count",               1000, 0, 35000);
    H["Theta_pi"]     = new TH1F("Theta_pi",
        "#pi^{-} polar angle;#theta (deg);Count",               100, 0, 10);
    // Lab frame: centred near x ~ -875 mm due to 25 mrad crossing angle
    H["XYatZDC_pi"]   = new TH2F("XYatZDC_pi",
        "Projected XY at ZDC (#pi^{-}) [lab];X (mm);Y (mm)",   100,-1200,-600, 100,-300, 300);
    // Beam-centred: beam axis subtracted, centred near (0,0)
    H["XYatZDC_pi_BC"] = new TH2F("XYatZDC_pi_BC",
        "Projected XY at ZDC (#pi^{-}) [beam-centred];X_{BC} (mm);Y_{BC} (mm)", 100,-300,300, 100,-300,300);
    H["ZvsTheta_pi"]  = new TH2F("ZvsTheta_pi",
        "Z from ZDC vs #theta (#pi^{-});#theta (deg);Z from ZDC (mm)", 100,0,10, 100,0,15000);
    H["ZvsRadius_pi"] = new TH2F("ZvsRadius_pi",
        "Z from ZDC vs R at ZDC (#pi^{-});R_{ZDC} (mm);Z from ZDC (mm)", 100,0,1500, 100,0,15000);

    // ── Daughter: neutron (PDG 2112) ────────────────────────────────────────
    H["Vertex_n"]   = new TH1F("Vertex_n",
        "Neutron vertex z;z_{vtx} (mm);Count",                 1000, 0, 35000);
    H["Endpoint_n"] = new TH1F("Endpoint_n",
        "Neutron endpoint z;z_{end} (mm);Count",               1000, 0, 35000);
    H["Theta_n"]    = new TH1F("Theta_n",
        "Neutron polar angle;#theta (deg);Count",               100, 0, 10);
    H["XYatZDC_n"]  = new TH2F("XYatZDC_n",
        "Projected XY at ZDC (n);X (mm);Y (mm)",               100,-1500,-200, 100,-500,500);

    // ── Daughter: pi0 (PDG 111) ─────────────────────────────────────────────
    H["Vertex_pi0"]   = new TH1F("Vertex_pi0",
        "#pi^{0} vertex z;z_{vtx} (mm);Count",                 1000, 0, 35000);
    H["Endpoint_pi0"] = new TH1F("Endpoint_pi0",
        "#pi^{0} endpoint z;z_{end} (mm);Count",               1000, 0, 35000);
    H["Theta_pi0"]    = new TH1F("Theta_pi0",
        "#pi^{0} polar angle;#theta (deg);Count",               100, 0, 10);
    H["XYatZDC_pi0"]  = new TH2F("XYatZDC_pi0",
        "Projected XY at ZDC (#pi^{0});X (mm);Y (mm)",         100,-1900,-400, 100,-800,800);

    // ── p + pi- correlations at the ZDC face ────────────────────────────────
    // Primary detector-design outputs:
    //
    //   DeltaR    — transverse distance between p and pi- at the detector face.
    //               The minimum of this distribution sets the required spatial
    //               resolution of the tracker.
    //
    //   DeltaRvsZ — DeltaR vs distance from the decay vertex to the ZDC.
    //               Lambdas decaying close to the ZDC produce small separations
    //               (hard to resolve); those decaying far upstream produce large
    //               separations (easy to resolve).
    //
    //   OpenAngle — opening angle between the p and pi- momentum directions.
    //               Together with DeltaRvsZ, constrains the angular resolution
    //               needed to reconstruct the decay vertex from two track
    //               directions measured at the detector plane.
    H["ThetaCorr"]    = new TH2F("ThetaCorr",
        "#theta p vs #pi^{-};#theta_{p} (deg);#theta_{#pi^{-}} (deg)", 100,0,5, 100,0,5);
    H["RadiusSumVsZ"] = new TH2F("RadiusSumVsZ",
        "R_{p}+R_{#pi^{-}} vs Z from ZDC;R_{p}+R_{#pi^{-}} (mm);Z from ZDC (mm)",
        100, 1000, 3000, 100, 0, 15000);
    H["DeltaR"]    = new TH1F("DeltaR",
        "Transverse separation p-#pi^{-} at ZDC;#DeltaR (mm);Count",    200, 0, 250);
    H["OpenAngle"] = new TH1F("OpenAngle",
        "Opening angle p-#pi^{-};#alpha (deg);Count",                    100, 0, 2);
    H["DeltaRvsZ"] = new TH2F("DeltaRvsZ",
        "#DeltaR vs Z from ZDC (key plot);Z from ZDC (mm);#DeltaR (mm)", 100,0,15000, 100,0,250);

    // ── DeltaR at each tracker plane (truth, scenario-independent) ───────────
    for (int ip = 1; ip <= Tracker::N_PLANES; ++ip) {
        string key = Form("DeltaR_pl%d", ip);
        double zpl  = plane_z(ip);
        H[key] = new TH1F(key.c_str(),
            Form("#DeltaR at plane %d (z=%.0f mm);#DeltaR (mm);Count", ip, zpl),
            200, 0, 250);
    }

    // ── DeltaR zoomed: region where tracker pitch matters (0 to 3 mm) ────────
    // The 0-250 mm histograms above have 1.25 mm/bin — the entire sub-mm region
    // is buried in the first bin. This set zooms into [0, 3 mm] with 0.05 mm/bin
    // (= finest pitch scenario) so the vertical sigma_x lines are clearly separated.
    // Overlaid in one canvas to show how the separation grows from plane 1 to plane 4.
    for (int ip = 1; ip <= Tracker::N_PLANES; ++ip) {
        string key = Form("DeltaR_pl%d_zoom", ip);
        double zpl  = plane_z(ip);
        H[key] = new TH1F(key.c_str(),
            Form("#DeltaR plane %d (z=%.0f mm) low-#DeltaR region"
                 ";#DeltaR (mm);Count", ip, zpl),
            60, 0, 3.0);
    }

    // ── Tracker resolution study — analytical ────────────────────────────────
    for (int i = 0; i < Tracker::N_SC; ++i) {
        string t = Form("_sc%d", i);
        H["SigmaZ"         +t] = new TH1F(("SigmaZ"         +t).c_str(),
            "#sigma_{z} (analytical);#sigma_{z} (mm);Events",
            200, 0, 5000);
        H["SigmaAlpha"     +t] = new TH1F(("SigmaAlpha"     +t).c_str(),
            "#sigma_{#alpha} (analytical);#sigma_{#alpha} (mrad);Events",
            100, 0, 5);
        H["SigmaZ_vsAlpha" +t] = new TH2F(("SigmaZ_vsAlpha" +t).c_str(),
            "#sigma_{z} vs #alpha_{true};#alpha_{true} (mrad);#sigma_{z} (mm)",
            30, 0, 15.0, 200, 0, 5000);
        H["SigmaZ_vsD"     +t] = new TH2F(("SigmaZ_vsD"     +t).c_str(),
            "#sigma_{z} vs D_{eff};D_{eff} (mm);#sigma_{z} (mm)",
            100, 0, 15000, 200, 0, 5000);
    }

    // ── Tracker resolution study — smearing ──────────────────────────────────
    for (int i = 0; i < Tracker::N_SC; ++i) {
        string t = Form("_sc%d", i);
        H["ResZ"           +t] = new TH1F(("ResZ"           +t).c_str(),
            "z_{reco}-z_{true};#Deltaz (mm);Events",
            400, -4000, 4000);
        H["ResAlpha"       +t] = new TH1F(("ResAlpha"       +t).c_str(),
            "#alpha_{reco}-#alpha_{true};#Delta#alpha (mrad);Events",
            200, -5, 5);
        H["ResZ_vsAlpha"   +t] = new TH2F(("ResZ_vsAlpha"   +t).c_str(),
            "#Deltaz vs #alpha_{true};#alpha_{true} (mrad);#Deltaz (mm)",
            30, 0, 15.0, 400, -4000, 4000);
        // Primary deliverable: reconstructed vertex z from smearing PCA.
        H["ZReco"          +t] = new TH1F(("ZReco"          +t).c_str(),
            "Reconstructed vertex z;z_{reco} (mm);Events",
            300, 15000, 40000);

        // Precision-based efficiency histograms.
        // These count events where |z_reco - z_true| < threshold,
        // directly comparable to the analytical efficiency (sigma_z < threshold).
        // One histogram per threshold — filled with 1 (pass) or 0 (fail)
        // so GetMean() gives the efficiency directly.
        for (int j = 0; j < Tracker::N_THRESH; ++j) {
            string key = Form("SmearEff_thr%d", j) + t;
            H[key] = new TH1F(key.c_str(),
                Form("|#Deltaz| < %.0f mm;pass (1) or fail (0);Events",
                     Tracker::EFF_THRESH[j]),
                2, -0.5, 1.5);
        }

        // ── Validation profile: RMS of Delta_z vs alpha ───────────────────
        // TProfile with option "S" stores the RMS (spread) of the y values
        // in each alpha bin, not just the mean. This allows direct comparison
        // with the analytical sigma_z curve from hSigmaZ_vsAlpha:
        //   - Profile RMS  = actual residual width from smearing (what the
        //                    reconstruction achieves in practice)
        //   - Analytical   = predicted width from sigma_z formula
        // Agreement validates the Gaussian approximation. Disagreement at
        // small alpha reveals where the Taylor expansion breaks down.
        H["ResZ_rms_vsAlpha" +t] = new TProfile(
            ("ResZ_rms_vsAlpha" +t).c_str(),
            "RMS(#Deltaz) vs #alpha_{true}"
            ";#alpha_{true} (mrad);RMS(#Deltaz) (mm)",
            30, 0, 15.0, "S");  // 30 bins x 0.5 mrad, range 0-15 mrad
    }
}

// ============================================================================
// Utility functions
// ============================================================================

// Polar angle between the momentum vector and the z-axis [degrees].
// Standard definition: theta = acos(pz / |p|).
double PolarAngle(double px, double py, double pz) {
    double p = std::sqrt(px*px + py*py + pz*pz);
    return std::acos(pz / p) * 180.0 / M_PI;
}

// Azimuthal angle in the transverse plane [degrees], range (-180, +180].
double PhiAngle(double px, double py) {
    return std::atan2(py, px) * 180.0 / M_PI;
}

// Opening angle between two momentum directions [degrees].
// Computed from the dot product of the two unit vectors in spherical
// coordinates. The clamp to [-1, 1] prevents NaN from acos() due to
// floating-point rounding pushing the dot product slightly out of range.
double OpeningAngle(const DaughterKin& a, const DaughterKin& b) {
    double dot =
        std::sin(a.theta*M_PI/180.0) * std::cos(a.phi*M_PI/180.0) *
        std::sin(b.theta*M_PI/180.0) * std::cos(b.phi*M_PI/180.0)
      + std::sin(a.theta*M_PI/180.0) * std::sin(a.phi*M_PI/180.0) *
        std::sin(b.theta*M_PI/180.0) * std::sin(b.phi*M_PI/180.0)
      + std::cos(a.theta*M_PI/180.0) * std::cos(b.theta*M_PI/180.0);
    dot = std::max(-1.0, std::min(1.0, dot));
    return std::acos(dot) * 180.0 / M_PI;
}

// ============================================================================
// process_event
//
// Called once per event. Loops over all MC particles, selects Lambdas
// according to the cuts in namespace Cuts, and fills histograms.
// ============================================================================
void process_event(const podio::Frame& event, int event_number) {
    try {
        const edm4hep::MCParticleCollection& mcParticles =
            event.get<edm4hep::MCParticleCollection>("MCParticles");

        for (int idx = 0; idx < (int)mcParticles.size(); ++idx) {
            edm4hep::MCParticle particle = mcParticles.at(idx);

            // ── All Lambdas with 2 parents ────────────────────────────────
            // 2-parent condition selects Lambdas from the primary interaction
            // vertex (parents = beam electron + proton), not from secondary
            // re-interactions inside the detector material.
            // NOTE: this block fills LambdaEnd histograms for ALL such Lambdas,
            // including secondary ones — this is intentional for diagnostics.
            // The primary Lambda (index=4) is selected later.
            if (particle.getPDG() == Cuts::LAMBDA_PDG &&
                particle.getParents().size() == 2)
            {
                double ez = particle.getEndpoint().z;
                h1("LambdaEnd")->Fill(ez);
                if      (ez <  6000)                               h1("LambdaEnd_1")->Fill(ez);
                else if (ez >= 6000  && ez <= 20800)               h1("LambdaEnd_2")->Fill(ez);
                // LambdaEnd_3: primary Lambda only (index=4), has daughters (not escaped),
                // decay in field-free window. This makes it directly comparable to
                // lambda_field_free and LambdaEnd_analysis.
                else if (ez >  20800 && ez <  Cuts::DECAY_Z_MAX
                         && particle.getDaughters().size() > 0
                         && particle.getObjectID().index == Cuts::LAMBDA_PRIMARY_IDX)
                    h1("LambdaEnd_3")->Fill(ez);

                edm4hep::Vector3d mom = particle.getMomentum();
                h1("LambdaAngleGen")->Fill(PolarAngle(mom.x, mom.y, mom.z));
            }

            // ── All Lambdas regardless of parentage ───────────────────────
            if (particle.getPDG() == Cuts::LAMBDA_PDG) {
                double vz = particle.getVertex().z;
                h1("LambdaGenZ")->Fill(vz);
                if      (vz <  6000)                h1("LambdaGen_1")->Fill(vz);
                else if (vz >= 6000 && vz <= 20800) h1("LambdaGen_2")->Fill(vz);
                else                                h1("LambdaGen_3")->Fill(vz);
            }

            // ── Primary Lambda only from here on ──────────────────────────
            // LAMBDA_PRIMARY_IDX (=4) is the object index the generator assigns
            // to the Lambda from the hard scatter. Secondary Lambdas from
            // re-interactions have different indices and are skipped.
            if (particle.getPDG()            != Cuts::LAMBDA_PDG ||
                particle.getObjectID().index != Cuts::LAMBDA_PRIMARY_IDX)
                continue;

            // ── Lambda with no daughters ───────────────────────────────────
            if (particle.getDaughters().size() == 0) {
                h1("LambdaZnoDaughters")->Fill(particle.getEndpoint().z);
                counters.lambda_escaped++;
                continue;
            }

            // ── Geometric decay vertex cut ─────────────────────────────────
            // Count each rejection reason separately so we can understand
            // the difference between LambdaEnd_3 and counters.lambda_charged.
            double ez = particle.getEndpoint().z;

            if (ez <= Cuts::DECAY_Z_MIN) {
                counters.lambda_before_magnet++;
                continue;
            }
            if (ez >= Cuts::DECAY_Z_MAX) {
                counters.lambda_beyond_zdc++;
                continue;
            }

            // Lambda is in the field-free window — count it regardless of channel
            counters.lambda_field_free++;

            // ── Channel selection ──────────────────────────────────────────
            // 2 daughters = charged channel (p + pi-) or neutral (n + pi0).
            // The PDG check in the daughter loop distinguishes them.
            // Count each channel separately before applying the cut.
            int ndau = (int)particle.getDaughters().size();
            if (ndau != Cuts::LAMBDA_N_DAUGHTERS) {
                counters.lambda_other_ndau++;
                continue;
            }

            // Check daughter PDGs to distinguish charged from neutral channel.
            // We do this before the main daughter loop so the counter is right.
            bool has_proton = false, has_pion_minus = false;
            bool has_neutron = false, has_pion_zero = false;
            for (const edm4hep::MCParticle& d : particle.getDaughters()) {
                int pdg = d.getPDG();
                if      (pdg ==  Cuts::PDG_PROTON)     has_proton     = true;
                else if (pdg ==  Cuts::PDG_PION_MINUS)  has_pion_minus = true;
                else if (pdg ==  Cuts::PDG_NEUTRON)     has_neutron    = true;
                else if (pdg ==  Cuts::PDG_PION_ZERO)   has_pion_zero  = true;
            }

            if (has_neutron || has_pion_zero) {
                // Neutral channel: n + pi0 — not used in tracking analysis
                counters.lambda_neutral++;
                // Still fill daughter histograms for the neutral channel
                double xv = particle.getEndpoint().x;
                double yv = particle.getEndpoint().y;
                double zv = ez;
                for (const edm4hep::MCParticle& daughter : particle.getDaughters()) {
                    int pdg = daughter.getPDG();
                    edm4hep::Vector3d mom = daughter.getMomentum();
                    if (pdg == Cuts::PDG_NEUTRON) {
                        DaughterKin kin_n;
                        kin_n.compute(mom.x, mom.y, mom.z, xv, yv, zv,
                                      Cuts::ZDC_Z, Cuts::CROSSING_ANGLE);
                        h1("Vertex_n")  ->Fill(daughter.getVertex().z);
                        h1("Endpoint_n")->Fill(daughter.getEndpoint().z);
                        h1("Theta_n")   ->Fill(kin_n.theta);
                        h2("XYatZDC_n") ->Fill(kin_n.xAtZDC, kin_n.yAtZDC);
                        counters.neutron++;
                    } else if (pdg == Cuts::PDG_PION_ZERO) {
                        DaughterKin kin_pi0;
                        kin_pi0.compute(mom.x, mom.y, mom.z, xv, yv, zv,
                                        Cuts::ZDC_Z, Cuts::CROSSING_ANGLE);
                        h1("Vertex_pi0")  ->Fill(daughter.getVertex().z);
                        h1("Endpoint_pi0")->Fill(daughter.getEndpoint().z);
                        h1("Theta_pi0")   ->Fill(kin_pi0.theta);
                        h2("XYatZDC_pi0") ->Fill(kin_pi0.xAtZDC, kin_pi0.yAtZDC);
                        counters.pion_zero++;
                    }
                }
                continue;
            }

            // Charged channel: p + pi- — used in the analysis
            counters.lambda_charged++;

            // Fill the analysis-sample endpoint histogram
            h1("LambdaEnd_analysis")->Fill(ez);

            // Lambda decay vertex position [mm].
            // All daughters originate from this point. The vertex is offset
            // from the z-axis by ~ez * tan(25 mrad) in x due to the beam
            // crossing angle — this must be included in the ZDC projection.
            double xv = particle.getEndpoint().x;
            double yv = particle.getEndpoint().y;
            double zv = particle.getEndpoint().z;   // == ez

            // Per-event kinematics for the proton and pi-.
            // Declared outside the daughter loop so they are accessible in
            // the correlation block below. isSet() returns false until
            // compute() is called for the matching PDG.
            DaughterKin kin_p, kin_pi;

            // ── Daughter loop ─────────────────────────────────────────────
            // Iterates over the MC particles produced at the Lambda decay
            // vertex. For the charged channel: proton (2212) + pi- (-211).
            for (const edm4hep::MCParticle& daughter : particle.getDaughters()) {
                int pdg = daughter.getPDG();
                edm4hep::Vector3d mom = daughter.getMomentum();

                if (pdg == Cuts::PDG_PROTON) {
                    kin_p.compute(mom.x, mom.y, mom.z, xv, yv, zv, Cuts::ZDC_Z, Cuts::CROSSING_ANGLE);
                    h1("Vertex_p")    ->Fill(daughter.getVertex().z);
                    h1("Endpoint_p")  ->Fill(daughter.getEndpoint().z);
                    h1("Theta_p")     ->Fill(kin_p.theta);
                    h2("XYatZDC_p")   ->Fill(kin_p.xAtZDC,  kin_p.yAtZDC);   // lab frame
                    h2("XYatZDC_p_BC")->Fill(kin_p.xBC,     kin_p.yBC);       // beam-centred
                    h2("ZvsTheta_p")  ->Fill(kin_p.theta,   Cuts::ZDC_Z - zv);
                    h2("ZvsRadius_p") ->Fill(kin_p.radius,  Cuts::ZDC_Z - zv);
                    counters.proton++;
                }
                else if (pdg == Cuts::PDG_PION_MINUS) {
                    kin_pi.compute(mom.x, mom.y, mom.z, xv, yv, zv, Cuts::ZDC_Z, Cuts::CROSSING_ANGLE);
                    h1("Vertex_pi")   ->Fill(daughter.getVertex().z);
                    h1("Endpoint_pi") ->Fill(daughter.getEndpoint().z);
                    h1("Theta_pi")    ->Fill(kin_pi.theta);
                    h2("XYatZDC_pi")   ->Fill(kin_pi.xAtZDC, kin_pi.yAtZDC);  // lab frame
                    h2("XYatZDC_pi_BC")->Fill(kin_pi.xBC,    kin_pi.yBC);      // beam-centred
                    h2("ZvsTheta_pi") ->Fill(kin_pi.theta,    Cuts::ZDC_Z - zv);
                    h2("ZvsRadius_pi")->Fill(kin_pi.radius,   Cuts::ZDC_Z - zv);
                    counters.pion_minus++;
                }
                // Note: neutron and pi0 daughters are handled in the neutral
                // channel block above and will not appear here because we only
                // reach this loop for charged channel Lambdas (has_proton=true).
            } // end daughter loop

            // ── p + pi- correlation histograms ────────────────────────────
            // Filled only when both daughters were found (isSet() == true).
            // See comments in book_histograms() for the physics motivation.
            if (kin_p.isSet() && kin_pi.isSet()) {
                double dist2ZDC = Cuts::ZDC_Z - zv;
                double deltaX = kin_p.xAtZDC - kin_pi.xAtZDC;
                double deltaY = kin_p.yAtZDC - kin_pi.yAtZDC;
                double deltaR = std::sqrt(deltaX*deltaX + deltaY*deltaY);
                double alpha  = OpeningAngle(kin_p, kin_pi);  // degrees

                h2("ThetaCorr")   ->Fill(kin_p.theta, kin_pi.theta);
                h2("RadiusSumVsZ")->Fill(kin_p.radius + kin_pi.radius, dist2ZDC);
                h1("DeltaR")      ->Fill(deltaR);
                h1("OpenAngle")   ->Fill(alpha);
                h2("DeltaRvsZ")   ->Fill(dist2ZDC, deltaR);

                // ── Tracker resolution study ──────────────────────────────
                // Physics goal: can we say "two tracks from a common vertex
                // in the field-free region" — consistent with Lambda decay?
                //
                // For each sigma_x scenario:
                //   1. Find first plane k where ΔR(z_k) > sigma_x
                //      ΔR scales linearly from vertex: ΔR(z) = ΔR_ZDC*(z-zv)/dist2ZDC
                //   2. D_eff = z_k - zv
                //   3. N_eff = N_PLANES - k + 1  (resolved planes available)
                //   4. L_eff = plane_z(N_PLANES) - plane_z(k)  (lever arm)
                //   5. sigma_theta_eff uses N_eff, L_eff — NOT the fixed N=4
                //      Key point: fewer resolved planes = worse angular resolution
                //      even if D_eff happens to be larger
                //   6. N_eff < 2: unresolvable (conservative — see note below)
                //
                // When all 4 planes resolve the tracks (first_plane=1):
                //   N_eff=4, L_eff=2000mm — best possible angular resolution
                //
                // Conservative assumption for N_eff=1: merged hits at earlier
                // planes constrain the Lambda direction but not individual track
                // slopes. A dedicated algorithm could recover some of these
                // events but is beyond the scope of this study.

                double alpha_rad = alpha * M_PI / 180.0;  // rad

                // ── ΔR at each plane: two methods compared ────────────────
                //
                // Method A (scaling): uses DeltaR already computed at ZDC face
                //   DeltaR(z_plane) = DeltaR_ZDC * (z_plane - zv) / dist2ZDC
                //   Derivation: both tracks are straight lines from (xv,yv,zv).
                //   The separation DeltaR grows linearly with (z - zv) because
                //   DeltaR(z) = ||(slope_p - slope_pi)|| * (z - zv) and the
                //   slope difference is constant. DeltaR_ZDC is at z=ZDC_Z so
                //   the ratio (z_plane-zv)/(ZDC_Z-zv) gives the scaling factor.
                //
                // Method B (direct): computes hit positions at the plane from
                //   the true track slopes tx=px/pz, ty=py/pz:
                //   x_p(z)  = xv + tx_p  * (z_plane - zv)
                //   x_pi(z) = xv + tx_pi * (z_plane - zv)
                //   DeltaR  = sqrt((x_p-x_pi)^2 + (y_p-y_pi)^2)
                //   This is the direct derivation — no reference to ZDC needed.
                //
                // Both methods are mathematically equivalent for straight-line
                // tracks. The residual (A - B) should be zero to floating-point
                // precision. Any non-zero residual would indicate a bug.

                // Compute slopes once for method B (reused in smearing below)
                double tx_p_true  = kin_p.theta  > 0 ?
                    std::tan(kin_p.theta *M_PI/180.0)*std::cos(kin_p.phi *M_PI/180.0) : 0.0;
                double ty_p_true  = kin_p.theta  > 0 ?
                    std::tan(kin_p.theta *M_PI/180.0)*std::sin(kin_p.phi *M_PI/180.0) : 0.0;
                double tx_pi_true = kin_pi.theta > 0 ?
                    std::tan(kin_pi.theta*M_PI/180.0)*std::cos(kin_pi.phi*M_PI/180.0) : 0.0;
                double ty_pi_true = kin_pi.theta > 0 ?
                    std::tan(kin_pi.theta*M_PI/180.0)*std::sin(kin_pi.phi*M_PI/180.0) : 0.0;

                for (int ip = 1; ip <= Tracker::N_PLANES; ++ip) {
                    double zpl = plane_z(ip);
                    if (zpl <= zv) continue;
                    double dR = deltaR * (zpl - zv) / dist2ZDC;
                    h1(Form("DeltaR_pl%d",      ip))->Fill(dR);
                    h1(Form("DeltaR_pl%d_zoom", ip))->Fill(dR);
                }

                if (alpha_rad > 1e-6) {
                    // Slopes tx_p_true etc. already computed above in the
                    // DeltaR comparison block — reused here directly.
                    for (int i = 0; i < Tracker::N_SC; ++i) {
                        string t = Form("_sc%d", i);

                        // Find first plane where tracks are resolved
                        int    first_plane = -1;
                        double D_eff       = -1.0;
                        for (int ip = 1; ip <= Tracker::N_PLANES; ++ip) {
                            double zpl = plane_z(ip);
                            if (zpl <= zv) continue;
                            if (deltaR * (zpl - zv) / dist2ZDC > Tracker::SIGMA_X[i]) {
                                first_plane = ip;
                                D_eff       = zpl - zv;
                                break;
                            }
                        }
                        if (first_plane < 0) continue; // no plane resolves tracks

                        // Angular resolution from resolved planes only
                        double s_th = sigma_theta_eff(Tracker::SIGMA_X[i], first_plane);
                        if (s_th > 1e8) continue; // N_eff < 2, unresolvable

                        // ── Analytical ────────────────────────────────────
                        double s_alpha = std::sqrt(2.0) * s_th;
                        double s_z     = std::sqrt(2.0) * D_eff * s_th / alpha_rad;
                        h1("SigmaZ"         +t)->Fill(s_z);
                        h1("SigmaAlpha"     +t)->Fill(s_alpha * 1e3);
                        h2("SigmaZ_vsAlpha" +t)->Fill(alpha_rad * 1e3, s_z);
                        h2("SigmaZ_vsD"     +t)->Fill(D_eff, s_z);

                        // ── Smearing: PCA from first resolved plane ───────
                        double z_ref    = plane_z(first_plane);
                        double dz_ref   = z_ref - zv;
                        double x_p_ref  = xv + dz_ref * tx_p_true;
                        double y_p_ref  = yv + dz_ref * ty_p_true;
                        double x_pi_ref = xv + dz_ref * tx_pi_true;
                        double y_pi_ref = yv + dz_ref * ty_pi_true;

                        double tx_p_reco  = tx_p_true  + gRandom->Gaus(0.0, s_th);
                        double ty_p_reco  = ty_p_true  + gRandom->Gaus(0.0, s_th);
                        double tx_pi_reco = tx_pi_true + gRandom->Gaus(0.0, s_th);
                        double ty_pi_reco = ty_pi_true + gRandom->Gaus(0.0, s_th);

                        double norm_p  = std::sqrt(tx_p_reco*tx_p_reco   + ty_p_reco*ty_p_reco   + 1.0);
                        double norm_pi = std::sqrt(tx_pi_reco*tx_pi_reco + ty_pi_reco*ty_pi_reco + 1.0);
                        double dot_reco =
                            (tx_p_reco*tx_pi_reco + ty_p_reco*ty_pi_reco + 1.0) / (norm_p * norm_pi);
                        dot_reco = std::max(-1.0, std::min(1.0, dot_reco));
                        double alpha_reco = std::acos(dot_reco);

                        double z_reco = zv; // fallback — not a real reconstruction
                        bool   reco_ok = false;
                        double dtx    = tx_p_reco - tx_pi_reco;
                        double dty    = ty_p_reco - ty_pi_reco;
                        if (std::abs(dtx) > 1e-9 && std::abs(dty) > 1e-9) {
                            z_reco   = 0.5 * (z_ref - (x_p_ref - x_pi_ref) / dtx
                                            + z_ref - (y_p_ref - y_pi_ref) / dty);
                            reco_ok  = true;
                        } else if (std::abs(dtx) > 1e-9) {
                            z_reco   = z_ref - (x_p_ref - x_pi_ref) / dtx;
                            reco_ok  = true;
                        } else if (std::abs(dty) > 1e-9) {
                            z_reco   = z_ref - (y_p_ref - y_pi_ref) / dty;
                            reco_ok  = true;
                        }

                        h1("ResZ"         +t)->Fill(z_reco - zv);
                        h1("ResAlpha"     +t)->Fill((alpha_reco - alpha_rad) * 1e3);
                        h2("ResZ_vsAlpha" +t)->Fill(alpha_rad * 1e3, z_reco - zv);
                        // Only fill the validation profile when the reconstruction
                        // actually produced a result — skip the fallback z_reco=zv
                        // case which would give a spurious zero residual.
                        if (reco_ok) {
                            static_cast<TProfile*>(H.at("ResZ_rms_vsAlpha" +t))
                                ->Fill(alpha_rad * 1e3, z_reco - zv);
                        }
                        // Primary deliverable
                        h1("ZReco" +t)->Fill(z_reco);

                        // Precision-based efficiency: |z_reco - z_true| < threshold
                        // Same thresholds as analytical study — directly comparable.
                        if (reco_ok) {
                            double res_z = std::abs(z_reco - zv);
                            for (int j = 0; j < Tracker::N_THRESH; ++j) {
                                string key = Form("SmearEff_thr%d", j) + t;
                                h1(key)->Fill(res_z < Tracker::EFF_THRESH[j] ? 1.0 : 0.0);
                            }
                        }
                    }
                }
                // ── End tracker study ─────────────────────────────────────
            }

        } // end particle loop

    } catch (const std::exception& e) {
        fmt::print("Error accessing MCParticles in event {}: {}\n",
                   event_number, e.what());
    }
}

// ============================================================================
// process_file
//
// Opens one PODIO ROOT file and calls process_event() for each entry, up to
// events_limit (all entries if events_limit < 0).
// ============================================================================
void process_file(const std::string& filename, int events_limit) {
    current_file = filename;
    try {
        podio::ROOTReader reader;
        reader.openFile(filename);

        unsigned int nEvents = reader.getEntries(podio::Category::Event);
        fmt::print("  {} events in file\n", nEvents);

        unsigned int limit = (events_limit < 0)
            ? nEvents
            : std::min((unsigned int)events_limit, nEvents);

        for (unsigned int i = 0; i < limit; ++i) {
            podio::Frame event(reader.readNextEntry(podio::Category::Event));
            process_event(event, (int)i);
            counters.total_events++;
        }
    } catch (const std::exception& e) {
        fmt::print("Error processing file {}: {}\n", filename, e.what());
    }
}

// ============================================================================
// plot
//
// Creates all display canvases and draws the histograms.
// Returns the canvas list so save_plots() can reuse them directly — this
// prevents the two functions from diverging as they did in the original code.
// ============================================================================
vector<TCanvas*> plot() {
    vector<TCanvas*> pages;

    // Helper: create a canvas, optionally divide it, push onto the page list.
    // Defined as a lambda to keep the canvas-creation lines below compact.
    auto make_canvas = [&](const char* name, const char* title,
                           int nx = 1, int ny = 1) -> TCanvas* {
        TCanvas* c = new TCanvas(name, title, 800, 600);
        if (nx > 1 || ny > 1) c->Divide(nx, ny);
        pages.push_back(c);
        return c;
    };

    // ── Lambda kinematics ─────────────────────────────────────────────────
    {
        TCanvas* c = make_canvas("cLambdaEnd", "Lambda decay endpoint", 2, 3);
        c->cd(1); h1("LambdaEnd")  ->Draw();
        c->cd(2); h1("LambdaEnd_1")->Draw();
        c->cd(3); h1("LambdaEnd_2")->Draw();
        c->cd(4); h1("LambdaEnd_3")->Draw();
        // Pad 5: overlay of all field-free Lambdas vs analysis sample.
        // Shows the neutral channel fraction as the gap between the two curves.
        c->cd(5);
        h1("LambdaEnd_3")->SetLineColor(kBlue); h1("LambdaEnd_3")->SetLineWidth(2);
        h1("LambdaEnd_3")->Draw("HIST");
        h1("LambdaEnd_analysis")->SetLineColor(kRed); h1("LambdaEnd_analysis")->SetLineWidth(2);
        h1("LambdaEnd_analysis")->Draw("HIST SAME");
        TLegend* leg = new TLegend(0.35, 0.65, 0.88, 0.85);
        leg->AddEntry(h1("LambdaEnd_3"),
            "All primary #Lambda, field-free (both channels)", "l");
        leg->AddEntry(h1("LambdaEnd_analysis"),
            "Analysis sample (p+#pi^{-} charged channel only)", "l");
        leg->Draw();
        // Pad 6: analysis sample alone, clean view of what feeds the tracker study
        c->cd(6); h1("LambdaEnd_analysis")->Draw("HIST");
        c->Update();
    }
    {
        TCanvas* c = make_canvas("cLambdaGen", "Lambda generation vertex", 2, 2);
        c->cd(1); h1("LambdaGenZ") ->Draw();
        c->cd(2); h1("LambdaGen_1")->Draw();
        c->cd(3); h1("LambdaGen_2")->Draw();
        c->cd(4); h1("LambdaGen_3")->Draw();
        c->Update();
    }
    {
        TCanvas* c = make_canvas("cLambdaMisc", "Lambda angle and no-daughters", 1, 2);
        c->cd(1); h1("LambdaAngleGen")    ->Draw();
        c->cd(2); h1("LambdaZnoDaughters")->Draw();
        c->Update();
    }

    // ── Proton ────────────────────────────────────────────────────────────
    {
        TCanvas* c = make_canvas("cProton", "Proton decay products", 2, 2);
        c->cd(1); h1("Theta_p")    ->Draw();
        c->cd(2); h2("XYatZDC_p") ->Draw("COLZ");  // lab frame
        c->cd(3); h2("ZvsRadius_p")->Draw("COLZ");
        c->cd(4); h2("ZvsTheta_p") ->Draw("COLZ");
        c->Update();
    }

    // ── Pi- ───────────────────────────────────────────────────────────────
    {
        TCanvas* c = make_canvas("cPiMinus", "#pi^{-} decay products", 2, 2);
        c->cd(1); h1("Theta_pi")    ->Draw();
        c->cd(2); h2("XYatZDC_pi") ->Draw("COLZ");  // lab frame
        c->cd(3); h2("ZvsRadius_pi")->Draw("COLZ");
        c->cd(4); h2("ZvsTheta_pi") ->Draw("COLZ");
        c->Update();
    }

    // ── Beam-centred XY: p and pi- relative to hadron beam axis ──────────
    // These show the spread seen by the detector independently of placement.
    // Both should be centred near (0,0) after subtracting the beam axis offset.
    {
        TCanvas* c = make_canvas("cXYbeamCentred",
            "XY at ZDC [beam-centred frame]", 2, 1);
        c->cd(1); h2("XYatZDC_p_BC") ->Draw("COLZ");
        c->cd(2); h2("XYatZDC_pi_BC")->Draw("COLZ");
        c->Update();
    }

    // ── Neutron ───────────────────────────────────────────────────────────
    {
        TCanvas* c = make_canvas("cNeutron", "Neutron decay products", 2, 2);
        c->cd(1); h1("Vertex_n")  ->Draw();
        c->cd(2); h1("Endpoint_n")->Draw();
        c->cd(3); h1("Theta_n")   ->Draw();
        c->cd(4); h2("XYatZDC_n") ->Draw("COLZ");
        c->Update();
    }

    // ── Pi0 ───────────────────────────────────────────────────────────────
    {
        TCanvas* c = make_canvas("cPiZero", "#pi^{0} decay products", 2, 2);
        c->cd(1); h1("Vertex_pi0")  ->Draw();
        c->cd(2); h1("Endpoint_pi0")->Draw();
        c->cd(3); h1("Theta_pi0")   ->Draw();
        c->cd(4); h2("XYatZDC_pi0") ->Draw("COLZ");
        c->Update();
    }

    // ── p + pi- correlations — key detector-design plots ─────────────────
    {
        TCanvas* c = make_canvas("cCorr", "p + #pi^{-} correlations", 3, 2);
        c->cd(1); h2("ThetaCorr")   ->Draw("COLZ");
        c->cd(2); h2("RadiusSumVsZ")->Draw("COLZ");
        c->cd(3); h1("DeltaR")      ->Draw();
        c->cd(4); h1("OpenAngle")   ->Draw();
        c->cd(5); h2("DeltaRvsZ")   ->Draw("COLZ");
        c->Update();
    }

    return pages;
}

// ============================================================================
// save_plots
//
// Writes all canvases from plot() into one multi-page PDF.
// ROOT multi-page PDF convention:
//   first page  -> canvas->Print("file.pdf(")
//   middle pages-> canvas->Print("file.pdf")
//   last page   -> canvas->Print("file.pdf)")
// ============================================================================
void save_plots(const vector<TCanvas*>& pages, const string& output_name) {
    if (pages.empty()) {
        fmt::print("save_plots: no pages to save.\n");
        return;
    }

    string pdf = output_name + ".pdf";
    fmt::print("Saving {} pages to {}\n", pages.size(), pdf);

    for (size_t i = 0; i < pages.size(); ++i) {
        if (!pages[i]) continue;

        bool first = (i == 0);
        bool last  = (i == pages.size() - 1);

        if      (first && last) {
            // Edge case: only one canvas total
            pages[i]->Print((pdf + "(").c_str());
            pages[i]->Print((pdf + ")").c_str());
        }
        else if (first) pages[i]->Print((pdf + "(").c_str());
        else if (last)  pages[i]->Print((pdf + ")").c_str());
        else            pages[i]->Print(pdf.c_str());
    }

    fmt::print("Done.\n");
}

// ============================================================================
// plot_tracker_study
//
// Draws analytical and smearing results in separate canvases.
// Returns pages for save_tracker_study.
// Efficiency calculation is embedded here: for each scenario and threshold,
// computes the fraction of events with sigma_z < threshold.
//
// Visual style: each scenario is drawn with a semi-transparent filled area
// (SetFillColorAlpha) plus a solid outline, making overlapping distributions
// clearly distinguishable without relying on colour alone.
// ============================================================================
vector<TCanvas*> plot_tracker_study() {
    vector<TCanvas*> pages;

    // Colours per scenario
    const int colors[Tracker::N_SC] = {
        kBlue, kRed, kGreen+2, kOrange+7, kViolet
    };

    // Fill styles: solid transparent fill with matching outline.
    // SetFillColorAlpha(color, alpha) — alpha=0.25 gives a light wash.
    // Draw option "HIST" for outline only; we call SetFillColorAlpha separately.
    auto style_hist = [&](TH1F* h, int isc) {
        h->SetLineColor(colors[isc]);
        h->SetLineWidth(2);
        h->SetFillColorAlpha(colors[isc], 0.25);
    };

    auto make_canvas = [&](const char* name, const char* title,
                           int nx = 1, int ny = 1) -> TCanvas* {
        TCanvas* c = new TCanvas(name, title, 800, 600);
        if (nx > 1 || ny > 1) c->Divide(nx, ny);
        pages.push_back(c);
        return c;
    };

    // Helper: draw one overlay canvas of 1D histograms from the map.
    // hname_prefix: e.g. "SigmaZ" — will look up "SigmaZ_sc0" .. "SigmaZ_sc4"
    // canvas_title: shown in the window title bar only (not the pad title)
    // pad_title: explicit title set on the pad after drawing
    // draw_logy: whether to set log y scale
    auto draw_overlay = [&](const char* hname_prefix,
                             const char* canvas_name,
                             const char* window_title,
                             const char* pad_title,
                             bool draw_logy) -> TCanvas* {
        TCanvas* c = make_canvas(canvas_name, window_title);
        if (draw_logy) c->SetLogy();
        TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string key = string(hname_prefix) + Form("_sc%d", i);
            TH1F* h = h1(key);
            style_hist(h, i);
            // "HIST" draws the outline; adding "SAME" for subsequent ones.
            // Use "HIST F" to also draw the fill.
            h->Draw(i == 0 ? "HIST F" : "HIST F SAME");
            leg->AddEntry(h, Tracker::LABELS[i], "f");
        }
        // Draw outlines on top so they are not hidden by fills
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string key = string(hname_prefix) + Form("_sc%d", i);
            h1(key)->Draw("HIST SAME");
        }
        leg->Draw();
        // Set the pad title explicitly — avoids ROOT inheriting the first
        // histogram's title (which includes the scenario label)
        c->GetPad(0)->SetTitle(pad_title);
        gPad->Modified();
        c->Update();
        return c;
    };

    // ── Page 1: DeltaR at each tracker plane (2x2) ───────────────────────
    {
        TCanvas* c = make_canvas("cDeltaR_planes",
            "#DeltaR at each tracker plane (truth)", 2, 2);
        const int pl_colors[4] = {kBlue, kRed, kGreen+2, kOrange+7};
        for (int ip = 1; ip <= Tracker::N_PLANES; ++ip) {
            c->cd(ip);
            TH1F* h = h1(Form("DeltaR_pl%d", ip));
            h->SetLineColor(pl_colors[ip-1]);
            h->SetLineWidth(2);
            h->Draw("HIST");
            for (int i = 0; i < Tracker::N_SC; ++i) {
                double sx = Tracker::SIGMA_X[i];
                if (sx > h->GetXaxis()->GetXmax()) continue;
                TLine* line = new TLine(sx, 0, sx, h->GetMaximum()*0.9);
                line->SetLineColor(colors[i]);
                line->SetLineStyle(2); line->SetLineWidth(1);
                line->Draw();
            }
        }
        c->Update();
    }

    // ── Page 2: DeltaR zoomed overlay — all planes, low-DeltaR region ────
    // Shows the fraction of events with DeltaR below each tracker pitch.
    // All four planes overlaid so the growth of separation across the
    // tracker is visible. Vertical lines mark the five sigma_x scenarios.
    // Events to the LEFT of a line are unresolved by that pitch.
    {
        TCanvas* c = make_canvas("cDeltaR_zoom",
            "#DeltaR zoomed: low-separation region vs tracker pitch");
        c->SetLogy();
        TLegend* legP = new TLegend(0.55, 0.55, 0.88, 0.88);
        legP->SetHeader("Tracker plane", "C");
        const int pl_colors[4] = {kBlue, kRed, kGreen+2, kOrange+7};
        for (int ip = 1; ip <= Tracker::N_PLANES; ++ip) {
            TH1F* h = h1(Form("DeltaR_pl%d_zoom", ip));
            h->SetLineColor(pl_colors[ip-1]);
            h->SetLineWidth(2);
            h->SetFillColorAlpha(pl_colors[ip-1], 0.15);
            h->Draw(ip == 1 ? "HIST F" : "HIST F SAME");
            legP->AddEntry(h, Form("Plane %d (z=%.0f mm)",
                ip, plane_z(ip)), "lf");
        }
        // Redraw outlines on top
        for (int ip = 1; ip <= Tracker::N_PLANES; ++ip)
            h1(Form("DeltaR_pl%d_zoom", ip))->Draw("HIST SAME");
        // Vertical lines at each sigma_x scenario with label
        TH1F* href = h1("DeltaR_pl1_zoom");
        double ymax = href->GetMaximum() * 2.0; // headroom for log scale
        TLegend* legS = new TLegend(0.12, 0.12, 0.42, 0.42);
        legS->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            double sx = Tracker::SIGMA_X[i];
            TLine* line = new TLine(sx, 0.5, sx, ymax);
            line->SetLineColor(colors[i]);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw();
            legS->AddEntry(line, Tracker::LABELS[i], "l");
        }
        legP->Draw();
        legS->Draw();
        c->Update();
    }

    // ── Analytical: page 1 — sigma_z overlay ─────────────────────────────
    draw_overlay("SigmaZ", "cSigmaZ",
        "Vertex z-resolution (analytical)",
        "Vertex z-resolution — all tracker pitches;#sigma_{z} (mm);Events",
        true);

    // ── Analytical: page 2 — sigma_alpha overlay ─────────────────────────
    draw_overlay("SigmaAlpha", "cSigmaAlpha",
        "Opening-angle resolution (analytical)",
        "Opening-angle resolution — all tracker pitches;#sigma_{#alpha} (mrad);Events",
        true);

    // ── Analytical: page 3 — sigma_z vs opening angle ────────────────────
    {
        TCanvas* c = make_canvas("cSigmaZ_vsAlpha",
            "#sigma_{z} vs opening angle (analytical)", 3, 2);
        for (int i = 0; i < Tracker::N_SC; ++i) {
            c->cd(i+1); gPad->SetLogz();
            h2("SigmaZ_vsAlpha" + string(Form("_sc%d",i)))->Draw("COLZ");
        }
        c->Update();
    }

    // ── Analytical: page 4 — sigma_z vs D ────────────────────────────────
    {
        TCanvas* c = make_canvas("cSigmaZ_vsD",
            "#sigma_{z} vs decay-to-tracker distance (analytical)", 3, 2);
        for (int i = 0; i < Tracker::N_SC; ++i) {
            c->cd(i+1); gPad->SetLogz();
            h2("SigmaZ_vsD" + string(Form("_sc%d",i)))->Draw("COLZ");
        }
        c->Update();
    }

    // ── Efficiency + summary: page 5 ─────────────────────────────────────
    {
        TCanvas* c = make_canvas("cSummary",
            "Mean #sigma_{z} and efficiency vs tracker pitch", 1, 2);

        std::vector<double> sx(Tracker::N_SC);
        std::vector<double> sy(Tracker::N_SC);
        std::vector<double> ey(Tracker::N_SC);
        double eff[Tracker::N_SC][Tracker::N_THRESH] = {};

        for (int i = 0; i < Tracker::N_SC; ++i) {
            string t   = Form("_sc%d", i);
            TH1F*  hSZ = h1("SigmaZ" + t);
            sx[i] = Tracker::SIGMA_X[i];
            sy[i] = hSZ->GetMean();
            ey[i] = hSZ->GetRMS();
            double total = hSZ->GetEntries();
            if (total > 0) {
                for (int j = 0; j < Tracker::N_THRESH; ++j) {
                    int bin = hSZ->FindBin(Tracker::EFF_THRESH[j]);
                    eff[i][j] = hSZ->Integral(1, bin) / total;
                }
            }
        }

        // Top pad: mean sigma_z vs sigma_x
        // Horizontal line at 1000 mm = primary threshold
        c->cd(1);
        gPad->SetLogx(); gPad->SetLogy();
        TGraphErrors* gMean = new TGraphErrors(
            Tracker::N_SC, sx.data(), sy.data(), nullptr, ey.data());
        gMean->SetTitle(
            "Mean #sigma_{z} vs tracker pitch"
            ";#sigma_{x} (mm);Mean #sigma_{z} (mm)");
        gMean->SetMarkerStyle(21);
        gMean->SetMarkerColor(kBlue);
        gMean->SetLineColor(kBlue);
        gMean->SetLineWidth(2);
        gMean->Draw("APL");
        // Primary threshold reference line at 1000 mm
        TLine* lThresh = new TLine(
            Tracker::SIGMA_X[0], Tracker::EFF_THRESH[Tracker::PRIMARY_THRESH],
            Tracker::SIGMA_X[Tracker::N_SC-1], Tracker::EFF_THRESH[Tracker::PRIMARY_THRESH]);
        lThresh->SetLineColor(kRed); lThresh->SetLineStyle(2); lThresh->SetLineWidth(2);
        lThresh->Draw();
        // Label
        TLatex* lbl = new TLatex(Tracker::SIGMA_X[0]*1.1,
            Tracker::EFF_THRESH[Tracker::PRIMARY_THRESH]*1.15,
            Form("%.0f mm target", Tracker::EFF_THRESH[Tracker::PRIMARY_THRESH]));
        lbl->SetTextColor(kRed); lbl->SetTextSize(0.04); lbl->Draw();

        // Bottom pad: efficiency curves per threshold
        // PRIMARY threshold (1000 mm) drawn thicker and labelled explicitly
        c->cd(2);
        gPad->SetLogx();
        // colors: 500mm=kBlue, 1000mm=kRed (primary), 2000mm=kGreen+2
        const int thr_colors[3] = {kBlue, kRed, kGreen+2};
        const int thr_width[3]  = {1, 3, 1}; // primary is thicker
        TLegend* legE = new TLegend(0.12, 0.12, 0.55, 0.42);
        legE->SetHeader("Vertex z threshold", "C");
        for (int j = 0; j < Tracker::N_THRESH; ++j) {
            std::vector<double> effv(Tracker::N_SC);
            for (int i = 0; i < Tracker::N_SC; ++i) effv[i] = eff[i][j];
            TGraph* gEff = new TGraph(Tracker::N_SC, sx.data(), effv.data());
            string label = Form("#sigma_{z} < %.0f mm", Tracker::EFF_THRESH[j]);
            if (j == Tracker::PRIMARY_THRESH) label += "  [PRIMARY]";
            gEff->SetTitle(
                "Efficiency vs tracker pitch"
                ";#sigma_{x} (mm);Fraction of events with #sigma_{z} < threshold");
            gEff->SetMarkerStyle(20 + j);
            gEff->SetMarkerSize(j == Tracker::PRIMARY_THRESH ? 1.2 : 0.8);
            gEff->SetMarkerColor(thr_colors[j]);
            gEff->SetLineColor(thr_colors[j]);
            gEff->SetLineWidth(thr_width[j]);
            gEff->Draw(j == 0 ? "APL" : "PL SAME");
            gEff->GetYaxis()->SetRangeUser(0.0, 1.05);
            legE->AddEntry(gEff, label.c_str(), "lp");
        }
        // 80% efficiency reference line — minimum acceptable performance
        TLine* l80 = new TLine(Tracker::SIGMA_X[0], 0.8,
                               Tracker::SIGMA_X[Tracker::N_SC-1], 0.8);
        l80->SetLineColor(kGray+2); l80->SetLineStyle(3); l80->SetLineWidth(1);
        l80->Draw();
        TLatex* lbl80 = new TLatex(Tracker::SIGMA_X[Tracker::N_SC-1]*0.7, 0.81, "80%");
        lbl80->SetTextColor(kGray+2); lbl80->SetTextSize(0.04); lbl80->Draw();
        legE->Draw();
        c->Update();

        // Print summary table — primary threshold highlighted
        fmt::print("\n========== Tracker Resolution Study ==========\n");
        fmt::print("  {} planes, L={:.0f} mm, standoff={:.0f} mm from ZDC\n",
                   Tracker::N_PLANES, Tracker::LENGTH, Tracker::STANDOFF);
        fmt::print("  PRIMARY threshold: {:.0f} mm (1 m vertex resolution)\n",
                   Tracker::EFF_THRESH[Tracker::PRIMARY_THRESH]);
        fmt::print("  {:>12s}  {:>16s}  {:>16s}  {:>18s}",
                   "sigma_x (mm)", "sigma_th (urad)", "sigma_al (mrad)",
                   "mean sigma_z (mm)");
        for (int j = 0; j < Tracker::N_THRESH; ++j)
            fmt::print("  {:>14s}", Form("eff<%.0fmm%s",
                Tracker::EFF_THRESH[j],
                j == Tracker::PRIMARY_THRESH ? "*" : ""));
        fmt::print("\n");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            double st = sigma_theta_eff(Tracker::SIGMA_X[i], 1) * 1e6;
            double sa = std::sqrt(2.0) * sigma_theta_eff(Tracker::SIGMA_X[i], 1) * 1e3;
            fmt::print("  {:>12.3f}  {:>16.1f}  {:>16.4f}  {:>14.1f}+/-{:.1f}",
                       Tracker::SIGMA_X[i], st, sa, sy[i], ey[i]);
            for (int j = 0; j < Tracker::N_THRESH; ++j)
                fmt::print("  {:>14.1f}%", eff[i][j] * 100.0);
            fmt::print("\n");
        }
        fmt::print("  (* = primary threshold)\n");
        fmt::print("==============================================\n\n");
    }

    // ── Page 8: ZReco distribution + field-free efficiency ───────────────
    // PRIMARY SMEARING DELIVERABLE.
    // Top: z_reco distributions for all scenarios. A good tracker places
    //      the reconstructed vertex inside [DECAY_Z_MIN, DECAY_Z_MAX].
    //      Tails outside = misreconstructed events.
    // Bottom: fraction of events with z_reco in the field-free window vs
    //         sigma_x — the honest efficiency from actual reconstruction.
    {
        TCanvas* c = make_canvas("cZReco",
            "Reconstructed vertex z — primary smearing deliverable", 1, 2);

        // Top pad: z_reco distributions overlaid
        c->cd(1);
        TLegend* legZ = new TLegend(0.12, 0.55, 0.45, 0.88);
        legZ->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string t = Form("_sc%d", i);
            TH1F* h = h1("ZReco" +t);
            style_hist(h, i);
            h->Draw(i == 0 ? "HIST F" : "HIST F SAME");
            legZ->AddEntry(h, Tracker::LABELS[i], "f");
        }
        for (int i = 0; i < Tracker::N_SC; ++i)
            h1("ZReco" + string(Form("_sc%d",i)))->Draw("HIST SAME");
        // Field-free region boundary lines
        double ymax = h1("ZReco_sc0")->GetMaximum();
        TLine* lmin = new TLine(Cuts::DECAY_Z_MIN, 0,
                                Cuts::DECAY_Z_MIN, ymax*0.95);
        TLine* lmax = new TLine(Cuts::DECAY_Z_MAX, 0,
                                Cuts::DECAY_Z_MAX, ymax*0.95);
        lmin->SetLineColor(kBlack); lmin->SetLineStyle(2); lmin->SetLineWidth(2);
        lmax->SetLineColor(kBlack); lmax->SetLineStyle(2); lmax->SetLineWidth(2);
        lmin->Draw(); lmax->Draw();
        legZ->Draw();

        // Bottom pad: precision-based smearing efficiency vs sigma_x
        // |z_reco - z_true| < threshold — directly comparable to analytical
        c->cd(2);
        gPad->SetLogx();
        const int thr_colors[3] = {kBlue, kRed, kGreen+2};
        const int thr_width[3]  = {1, 3, 1};
        TLegend* legE = new TLegend(0.12, 0.12, 0.55, 0.45);
        legE->SetHeader("|z_{reco}-z_{true}| threshold", "C");

        // Print comparison table: analytical vs smearing
        fmt::print("\n===== Efficiency comparison: analytical vs smearing =====\n");
        fmt::print("  (fraction of events with sigma_z or |Delta_z| < threshold)\n");
        fmt::print("  {:>12s}", "sigma_x");
        for (int j = 0; j < Tracker::N_THRESH; ++j)
            fmt::print("  {:>12s}  {:>12s}",
                Form("anal<%.0f", Tracker::EFF_THRESH[j]),
                Form("smear<%.0f", Tracker::EFF_THRESH[j]));
        fmt::print("\n");

        for (int j = 0; j < Tracker::N_THRESH; ++j) {
            std::vector<double> sx_s(Tracker::N_SC), eff_s(Tracker::N_SC);
            for (int i = 0; i < Tracker::N_SC; ++i) {
                string t   = Form("_sc%d", i);
                string key = Form("SmearEff_thr%d", j) + t;
                TH1F*  h   = h1(key);
                sx_s[i]    = Tracker::SIGMA_X[i];
                // bin 2 = "pass" (value=1), GetMean() gives fraction directly
                double total = h->GetEntries();
                eff_s[i] = total > 0 ?
                    h->GetBinContent(2) / total : 0.0;
            }
            TGraph* gEff = new TGraph(Tracker::N_SC, sx_s.data(), eff_s.data());
            string label = Form("|#Deltaz| < %.0f mm", Tracker::EFF_THRESH[j]);
            if (j == Tracker::PRIMARY_THRESH) label += "  [PRIMARY]";
            gEff->SetTitle(
                "Precision efficiency (smearing) vs tracker pitch"
                ";#sigma_{x} (mm);Fraction with |z_{reco}-z_{true}| < threshold");
            gEff->SetMarkerStyle(20 + j);
            gEff->SetMarkerSize(j == Tracker::PRIMARY_THRESH ? 1.2 : 0.8);
            gEff->SetMarkerColor(thr_colors[j]);
            gEff->SetLineColor(thr_colors[j]);
            gEff->SetLineWidth(thr_width[j]);
            gEff->Draw(j == 0 ? "APL" : "PL SAME");
            gEff->GetYaxis()->SetRangeUser(0.0, 1.05);
            legE->AddEntry(gEff, label.c_str(), "lp");
        }
        TLine* l80e = new TLine(Tracker::SIGMA_X[0], 0.8,
                                Tracker::SIGMA_X[Tracker::N_SC-1], 0.8);
        l80e->SetLineColor(kGray+2); l80e->SetLineStyle(3); l80e->Draw();
        legE->Draw();
        c->Update();

        // Print comparison table to stdout
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string t = Form("_sc%d", i);
            fmt::print("  {:>12s}", Tracker::LABELS[i]);
            for (int j = 0; j < Tracker::N_THRESH; ++j) {
                // Analytical
                TH1F* hSZ    = h1("SigmaZ" + t);
                int   bin    = hSZ->FindBin(Tracker::EFF_THRESH[j]);
                double anal  = hSZ->GetEntries() > 0 ?
                    hSZ->Integral(1, bin) / hSZ->GetEntries() : 0.0;
                // Smearing
                string key   = Form("SmearEff_thr%d", j) + t;
                TH1F*  hSE   = h1(key);
                double total = hSE->GetEntries();
                double smear = total > 0 ?
                    hSE->GetBinContent(2) / total : 0.0;
                fmt::print("  {:>11.1f}%  {:>11.1f}%", anal*100, smear*100);
            }
            fmt::print("\n");
        }
        fmt::print("=========================================================\n\n");

        c->Update();
    }

    // ── Smearing: page 9 — vertex z residuals overlay ────────────────────
    draw_overlay("ResZ", "cResZ",
        "Vertex z residual (smearing, tx/ty slopes)",
        "Vertex z residual — all tracker pitches;#Deltaz = z_{reco}-z_{true} (mm);Events",
        true);

    // ── Smearing: page 10 — opening angle residuals overlay ──────────────
    draw_overlay("ResAlpha", "cResAlpha",
        "Opening-angle residual (smearing, tx/ty slopes)",
        "Opening-angle residual — all tracker pitches;#Delta#alpha = #alpha_{reco}-#alpha_{true} (mrad);Events",
        true);

    // ── Smearing: page 11 — residual vs opening angle, per scenario ──────
    {
        TCanvas* c = make_canvas("cResZ_vsAlpha",
            "#Deltaz vs opening angle (smearing)", 3, 2);
        for (int i = 0; i < Tracker::N_SC; ++i) {
            c->cd(i+1);
            h2("ResZ_vsAlpha" + string(Form("_sc%d",i)))->Draw("COLZ");
        }
        c->Update();
    }

    // ── Validation: page 12 — RMS(Delta_z) vs alpha overlaid with sigma_z ─
    // Red points: mean analytical sigma_z per alpha bin from ProfileX of
    //             SigmaZ_vsAlpha. Each point = mean of sigma_z values in
    //             that alpha slice, error = error on the mean.
    // Blue points: RMS of Delta_z per alpha bin from ProjectionY of
    //             ResZ_vsAlpha. Each bin is projected onto the Delta_z axis,
    //             GetRMS() gives the standard deviation of the distribution
    //             in that alpha slice — directly comparable to sigma_z.
    //             This is the same quantity you see when you click a bin in
    //             the 2D histogram and look at the projection's Std Dev.
    {
        TCanvas* c = make_canvas("cValidation",
            "Validation: RMS(#Deltaz) vs analytical #sigma_{z}", 3, 2);

        for (int i = 0; i < Tracker::N_SC; ++i) {
            c->cd(i+1);
            string t = Form("_sc%d", i);

            // ── Red: analytical sigma_z mean per alpha bin ────────────────
            TProfile* profSigmaZ = h2("SigmaZ_vsAlpha"+t)->ProfileX(
                Form("profSigmaZ_sc%d", i));
            profSigmaZ->SetLineColor(kRed);
            profSigmaZ->SetLineWidth(2);
            profSigmaZ->SetMarkerColor(kRed);
            profSigmaZ->SetMarkerStyle(20);
            profSigmaZ->SetMarkerSize(0.8);
            profSigmaZ->SetTitle(
                Form("#sigma_{x}=%s;#alpha_{true} (mrad);"
                     "#sigma_{z} or RMS(#Deltaz) (mm)",
                     Tracker::LABELS[i]));
            profSigmaZ->Draw("PE");

            // ── Blue: RMS of Delta_z per alpha bin from ProjectionY ───────
            // Loop over x bins of ResZ_vsAlpha. For each bin project onto
            // the Delta_z axis and take GetRMS() of that 1D distribution.
            // This is identical to what you see when you manually project
            // a bin in the ROOT GUI — no approximation, no TProfile issues.
            TH2F* h2RZ = h2("ResZ_vsAlpha" +t);
            int   nbinsX = h2RZ->GetNbinsX();
            vector<double> gx, gy, gex, gey;
            for (int b = 1; b <= nbinsX; ++b) {
                TH1D* proj = h2RZ->ProjectionY(
                    Form("proj_sc%d_bin%d", i, b), b, b);
                double n = proj->GetEntries();
                if (n < 2) { delete proj; continue; }
                double rms   = proj->GetRMS();       // std dev of Delta_z
                double alpha = h2RZ->GetXaxis()->GetBinCenter(b);
                double width = h2RZ->GetXaxis()->GetBinWidth(b) * 0.5;
                // Uncertainty on RMS: sigma_RMS = RMS / sqrt(2*(n-1))
                double e_rms = rms / std::sqrt(2.0 * (n - 1.0));
                gx.push_back(alpha);
                gy.push_back(rms);
                gex.push_back(width);
                gey.push_back(e_rms);
                delete proj;
            }

            if (!gx.empty()) {
                TGraphErrors* grRMS_smear = new TGraphErrors(
                    (int)gx.size(),
                    gx.data(), gy.data(), gex.data(), gey.data());
                grRMS_smear->SetLineColor(kBlue);
                grRMS_smear->SetLineWidth(2);
                grRMS_smear->SetMarkerColor(kBlue);
                grRMS_smear->SetMarkerStyle(21);
                grRMS_smear->SetMarkerSize(0.8);
                grRMS_smear->Draw("PE SAME");

                if (i == 0) {
                    TLegend* leg = new TLegend(0.35, 0.65, 0.88, 0.88);
                    leg->AddEntry(profSigmaZ,
                        "Analytical #sigma_{z} (formula)", "lp");
                    leg->AddEntry(grRMS_smear,
                        "RMS(#Deltaz) from smearing", "lp");
                    leg->Draw();
                }
            }
        }
        c->Update();
    }

    // ── Page 13: Residual contamination after tracker veto ────────────────
    // Physics context:
    //   Without tracker: charged-channel p+pi- decays reaching the ZDC
    //   cannot be distinguished from neutral-channel n+pi0 decays.
    //   They contaminate the neutral-channel measurement.
    //
    //   Field-free sample (only clean Lambda decays, excluding material
    //   interactions):
    //     N_charged = 1324  (p + pi-, used in tracker study)
    //     N_neutral =  721  (n + pi0)
    //     N_other   =  166  (material interactions — excluded, not real decays)
    //     Total clean = 2045
    //
    //   Starting contamination without tracker:
    //     C0 = 1324 / 2045 = 64.7%
    //
    //   After tracker identifies fraction epsilon of charged events:
    //     C_residual = (1 - epsilon) * C0
    //
    //   Shown for both analytical (sigma_z < threshold) and smearing
    //   (|z_reco - z_true| < threshold) efficiency estimates.
    {
        const double N_charged = 1324.0;
        const double N_neutral =  721.0;
        const double C0 = N_charged / (N_charged + N_neutral);

        TCanvas* c = make_canvas("cContamination",
            "Residual charged-channel contamination after tracker veto", 2, 2);

        const int thr_colors[3] = {kBlue, kRed, kGreen+2};
        const int thr_width[3]  = {1, 3, 1};

        // Helper lambda to draw reference line and label
        auto draw_ref = [&](double yval) {
            TLine* l = new TLine(Tracker::SIGMA_X[0], yval,
                                 Tracker::SIGMA_X[Tracker::N_SC-1], yval);
            l->SetLineColor(kBlack); l->SetLineStyle(2); l->SetLineWidth(2);
            l->Draw();
            TLatex* lbl = new TLatex(Tracker::SIGMA_X[0]*1.1, yval*1.04,
                Form("No tracker: %.1f%%", yval));
            lbl->SetTextSize(0.04); lbl->Draw();
        };

        // ── Top left: analytical ──────────────────────────────────────────
        c->cd(1);
        gPad->SetLogx();
        {
            TLegend* leg = new TLegend(0.30, 0.50, 0.88, 0.88);
            leg->SetHeader("Analytical (#sigma_{z} < #tau)", "C");
            bool first = true;
            for (int j = 0; j < Tracker::N_THRESH; ++j) {
                std::vector<double> sx(Tracker::N_SC), cont(Tracker::N_SC);
                for (int i = 0; i < Tracker::N_SC; ++i) {
                    TH1F* h  = h1("SigmaZ" + string(Form("_sc%d",i)));
                    double n = h->GetEntries();
                    double e = n > 0 ?
                        h->Integral(1, h->FindBin(Tracker::EFF_THRESH[j]))/n : 0.0;
                    sx[i]   = Tracker::SIGMA_X[i];
                    cont[i] = (1.0 - e) * C0 * 100.0;
                }
                TGraph* g = new TGraph(Tracker::N_SC, sx.data(), cont.data());
                string lbl = Form("#tau=%.0f mm%s", Tracker::EFF_THRESH[j],
                    j==Tracker::PRIMARY_THRESH ? " [PRIMARY]" : "");
                g->SetTitle("Analytical;#sigma_{x} (mm);Residual contamination (%)");
                g->SetMarkerStyle(20+j); g->SetMarkerSize(0.9);
                g->SetMarkerColor(thr_colors[j]); g->SetLineColor(thr_colors[j]);
                g->SetLineWidth(thr_width[j]);
                g->Draw(first ? "APL" : "PL SAME");
                if (first) { g->GetYaxis()->SetRangeUser(0.0, C0*100.0*1.15); first=false; }
                leg->AddEntry(g, lbl.c_str(), "lp");
            }
            draw_ref(C0 * 100.0);
            leg->Draw();
        }

        // ── Top right: smearing ───────────────────────────────────────────
        c->cd(2);
        gPad->SetLogx();
        {
            TLegend* leg = new TLegend(0.30, 0.50, 0.88, 0.88);
            leg->SetHeader("Smearing (|#Deltaz| < #tau)", "C");
            bool first = true;
            for (int j = 0; j < Tracker::N_THRESH; ++j) {
                std::vector<double> sx(Tracker::N_SC), cont(Tracker::N_SC);
                for (int i = 0; i < Tracker::N_SC; ++i) {
                    string key  = Form("SmearEff_thr%d_sc%d", j, i);
                    TH1F*  h    = h1(key);
                    double n    = h->GetEntries();
                    double e    = n > 0 ? h->GetBinContent(2)/n : 0.0;
                    sx[i]   = Tracker::SIGMA_X[i];
                    cont[i] = (1.0 - e) * C0 * 100.0;
                }
                TGraph* g = new TGraph(Tracker::N_SC, sx.data(), cont.data());
                string lbl = Form("#tau=%.0f mm%s", Tracker::EFF_THRESH[j],
                    j==Tracker::PRIMARY_THRESH ? " [PRIMARY]" : "");
                g->SetTitle("Smearing;#sigma_{x} (mm);Residual contamination (%)");
                g->SetMarkerStyle(20+j); g->SetMarkerSize(0.9);
                g->SetMarkerColor(thr_colors[j]); g->SetLineColor(thr_colors[j]);
                g->SetLineWidth(thr_width[j]);
                g->Draw(first ? "APL" : "PL SAME");
                if (first) { g->GetYaxis()->SetRangeUser(0.0, C0*100.0*1.15); first=false; }
                leg->AddEntry(g, lbl.c_str(), "lp");
            }
            draw_ref(C0 * 100.0);
            leg->Draw();
        }

        // ── Bottom left: primary threshold, analytical vs smearing ────────
        c->cd(3);
        gPad->SetLogx();
        {
            int j = Tracker::PRIMARY_THRESH;
            std::vector<double> sx(Tracker::N_SC);
            std::vector<double> cA(Tracker::N_SC), cS(Tracker::N_SC);
            for (int i = 0; i < Tracker::N_SC; ++i) {
                sx[i] = Tracker::SIGMA_X[i];
                TH1F* hA = h1("SigmaZ" + string(Form("_sc%d",i)));
                double nA = hA->GetEntries();
                double eA = nA > 0 ?
                    hA->Integral(1, hA->FindBin(Tracker::EFF_THRESH[j]))/nA : 0.0;
                cA[i] = (1.0 - eA) * C0 * 100.0;
                string key = Form("SmearEff_thr%d_sc%d", j, i);
                TH1F* hS   = h1(key);
                double nS  = hS->GetEntries();
                double eS  = nS > 0 ? hS->GetBinContent(2)/nS : 0.0;
                cS[i] = (1.0 - eS) * C0 * 100.0;
            }
            TGraph* gA = new TGraph(Tracker::N_SC, sx.data(), cA.data());
            TGraph* gS = new TGraph(Tracker::N_SC, sx.data(), cS.data());
            gA->SetTitle(Form("Primary #tau=%.0f mm: analytical vs smearing"
                ";#sigma_{x} (mm);Residual contamination (%%)",
                Tracker::EFF_THRESH[j]));
            gA->SetMarkerStyle(20); gA->SetMarkerColor(kRed);
            gA->SetLineColor(kRed); gA->SetLineWidth(3);
            gS->SetMarkerStyle(21); gS->SetMarkerColor(kBlue);
            gS->SetLineColor(kBlue); gS->SetLineWidth(3);
            gA->Draw("APL");
            gA->GetYaxis()->SetRangeUser(0.0, C0*100.0*1.15);
            gS->Draw("PL SAME");
            draw_ref(C0 * 100.0);
            TLegend* leg = new TLegend(0.35, 0.65, 0.88, 0.88);
            leg->AddEntry(gA, "Analytical", "lp");
            leg->AddEntry(gS, "Smearing",   "lp");
            leg->Draw();
        }

        // ── Bottom right: text summary ────────────────────────────────────
        c->cd(4);
        gPad->Clear();
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.055);
        tx.SetTextFont(42);
        double y = 0.92;
        auto line = [&](const char* s) { tx.DrawLatex(0.05, y, s); y -= 0.08; };
        line("Field-free sample (clean decays only):");
        line(Form("  N_{{charged}} = %.0f  (p+#pi^{{-}})", N_charged));
        line(Form("  N_{{neutral}} = %.0f  (n+#pi^{{0}})", N_neutral));
        line(Form("  N_{{other}}   = 166  (material, excluded)"));
        line(Form("  C_{{0}} (no tracker) = %.1f%%", C0*100.0));
        line("C_{{residual}} = (1-#varepsilon) #times C_{{0}}");
        c->Update();

        // Print contamination table to stdout
        fmt::print("\n===== Residual contamination after tracker veto =====\n");
        fmt::print("  N_charged={:.0f}  N_neutral={:.0f}  C0={:.1f}%\n",
                   N_charged, N_neutral, C0*100.0);
        fmt::print("  {:>10s}", "sigma_x");
        for (int j = 0; j < Tracker::N_THRESH; ++j)
            fmt::print("  {:>10s}  {:>10s}",
                Form("anal%.0f%%", Tracker::EFF_THRESH[j]),
                Form("smear%.0f%%", Tracker::EFF_THRESH[j]));
        fmt::print("\n");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            fmt::print("  {:>10s}", Tracker::LABELS[i]);
            for (int j = 0; j < Tracker::N_THRESH; ++j) {
                TH1F* hA  = h1("SigmaZ" + string(Form("_sc%d",i)));
                double nA = hA->GetEntries();
                double eA = nA > 0 ?
                    hA->Integral(1,hA->FindBin(Tracker::EFF_THRESH[j]))/nA : 0.0;
                string key = Form("SmearEff_thr%d_sc%d",j,i);
                TH1F*  hS  = h1(key);
                double nS  = hS->GetEntries();
                double eS  = nS > 0 ? hS->GetBinContent(2)/nS : 0.0;
                fmt::print("  {:>9.1f}%  {:>9.1f}%",
                    (1.0-eA)*C0*100.0, (1.0-eS)*C0*100.0);
            }
            fmt::print("\n");
        }
        fmt::print("=====================================================\n\n");
    }

    return pages;
}

// ============================================================================
// save_tracker_study — writes all tracker study pages to a dedicated PDF
// ============================================================================
void save_tracker_study(const vector<TCanvas*>& pages) {
    if (pages.empty()) return;
    string pdf = string(Tracker::OUTPUT_PDF) + ".pdf";
    fmt::print("Saving tracker study: {} pages to {}\n", pages.size(), pdf);
    for (size_t i = 0; i < pages.size(); ++i) {
        if (!pages[i]) continue;
        bool first = (i == 0), last = (i == pages.size()-1);
        if      (first && last) {
            pages[i]->Print((pdf+"(").c_str());
            pages[i]->Print((pdf+")").c_str());
        }
        else if (first) pages[i]->Print((pdf+"(").c_str());
        else if (last)  pages[i]->Print((pdf+")").c_str());
        else            pages[i]->Print(pdf.c_str());
    }
    fmt::print("Done.\n");
}

// ============================================================================
// save_histograms
//
// Writes every histogram in the global map H to a ROOT file.
// This allows reopening later to replot, refit, or overlay with other
// datasets without reprocessing the full event loop.
//
// Usage after running:
//   TFile* f = TFile::Open("lambda_analysis_18x275.root");
//   TH1F* h = (TH1F*)f->Get("DeltaR");
//   h->Draw();
// ============================================================================
void save_histograms(const string& rootfile) {
    TFile* f = TFile::Open(rootfile.c_str(), "RECREATE");
    if (!f || f->IsZombie()) {
        fmt::print("ERROR: could not open {} for writing\n", rootfile);
        return;
    }
    int n = 0;
    for (auto& kv : H) {
        if (kv.second) { kv.second->Write(); ++n; }
    }
    f->Close();
    fmt::print("Saved {} histograms to {}\n", n, rootfile);
}

// ============================================================================
// lambda_analysis — ROOT macro entry point
//
// events_limit < 0  : process all events in each file (default)
// events_limit >= 0 : process at most N events per file (for quick tests)
//
// Examples:
//   root [0] lambda_analysis()       // all events in all files
//   root [0] lambda_analysis(500)    // 500 events per file, quick check
// ============================================================================
void lambda_analysis(int events_limit = -1) {
    book_histograms();

    for (int i = 1; i <= Cuts::N_FILES; ++i) {
        string filename = fmt::format(Cuts::FILE_PATTERN, i);
        fmt::print("Processing file {:2d}/{}: {}\n", i, Cuts::N_FILES, filename);
        process_file(filename, events_limit);
    }

    counters.print();

    // ── Save all histograms to ROOT file ──────────────────────────────────
    // Allows reopening later to replot without rerunning the event loop.
    save_histograms("lambda_analysis_18x275.root");

    // ── Produce output PDFs ───────────────────────────────────────────────
    vector<TCanvas*> pages = plot();
    save_plots(pages, Cuts::OUTPUT_PDF);

    vector<TCanvas*> tracker_pages = plot_tracker_study();
    save_tracker_study(tracker_pages);

    fmt::print("\nOutputs written:\n");
    fmt::print("  lambda_analysis_18x275.root  (all histograms)\n");
    fmt::print("  {}.pdf  (diagnostic + particle plots)\n", Cuts::OUTPUT_PDF);
    fmt::print("  {}.pdf  (tracker resolution study)\n", Tracker::OUTPUT_PDF);
}



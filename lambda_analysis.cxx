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
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
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

    // Efficiency thresholds: fraction of events with sigma_z below these [mm]
    constexpr double EFF_THRESH[] = {500.0, 1000.0, 2000.0};
    constexpr int    N_THRESH     = 3;

    constexpr const char* OUTPUT_PDF = "tracker_study_18x275";
}

// Angular resolution [rad] for a straight-line fit through N equally-spaced
// planes over total length L with per-plane spatial resolution sigma_x [mm].
//   sigma_theta = (sigma_x / L) * sqrt(12*(N-1) / (N*(N+1)))
// For N=4: factor = sqrt(12*3/(4*5)) = sqrt(1.8) = 1.34164...
inline double sigma_theta_rad(double sigma_x_mm) {
    constexpr double factor = 1.3416407864998738; // sqrt(1.8), exact for N=4
    return (sigma_x_mm / Tracker::LENGTH) * factor;
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
    int total_events = 0;
    int lambda       = 0;   // Lambdas decaying in the acceptance window
    int proton       = 0;   // protons from those Lambdas
    int pion_minus   = 0;   // pi- from those Lambdas
    int neutron      = 0;   // neutrons (neutral decay channel)
    int pion_zero    = 0;   // pi0     (neutral decay channel)

    void print() const {
        fmt::print("\n========== Summary ==========\n");
        fmt::print("Total events processed : {}\n", total_events);
        fmt::print("Lambda (after magnet)  : {}\n", lambda);
        fmt::print("  -> proton            : {}\n", proton);
        fmt::print("  -> pi-               : {}\n", pion_minus);
        fmt::print("  -> neutron           : {}\n", neutron);
        fmt::print("  -> pi0               : {}\n", pion_zero);
        fmt::print("=============================\n");
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
        "Lambda decay endpoint (>20.8 m);z_{end} (mm);Count",  1000, 20800, 35000);

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

    // ── Tracker resolution study — analytical ────────────────────────────────
    // Per scenario: sigma_z distribution, sigma_alpha, and 2D correlations.
    for (int i = 0; i < Tracker::N_SC; ++i) {
        string t = Form("_sc%d", i);
        string l = Tracker::LABELS[i];
        H["SigmaZ"         +t] = new TH1F(("SigmaZ"         +t).c_str(),
            ("#sigma_{z} [#sigma_{x}="+l+"];#sigma_{z} (mm);Events").c_str(),
            200, 0, 5000);
        H["SigmaAlpha"     +t] = new TH1F(("SigmaAlpha"     +t).c_str(),
            ("#sigma_{#alpha} [#sigma_{x}="+l+"];#sigma_{#alpha} (mrad);Events").c_str(),
            100, 0, 5);
        H["SigmaZ_vsAlpha" +t] = new TH2F(("SigmaZ_vsAlpha" +t).c_str(),
            ("#sigma_{z} vs #alpha [#sigma_{x}="+l+
             "];#alpha_{true} (mrad);#sigma_{z} (mm)").c_str(),
            100, 0, 100, 200, 0, 5000);
        H["SigmaZ_vsD"     +t] = new TH2F(("SigmaZ_vsD"     +t).c_str(),
            ("#sigma_{z} vs D [#sigma_{x}="+l+
             "];D (mm);#sigma_{z} (mm)").c_str(),
            100, 0, 15000, 200, 0, 5000);
    }

    // ── Tracker resolution study — smearing ──────────────────────────────────
    // Residuals from reconstructed vs true quantities after Gaussian smearing
    // of the track angles. One set per scenario.
    for (int i = 0; i < Tracker::N_SC; ++i) {
        string t = Form("_sc%d", i);
        string l = Tracker::LABELS[i];
        // Vertex z residual: reco - true [mm]
        H["ResZ"           +t] = new TH1F(("ResZ"           +t).c_str(),
            ("z_{reco}-z_{true} [#sigma_{x}="+l+
             "];#Deltaz (mm);Events").c_str(),
            400, -4000, 4000);
        // Opening angle residual: reco - true [mrad]
        H["ResAlpha"       +t] = new TH1F(("ResAlpha"       +t).c_str(),
            ("#alpha_{reco}-#alpha_{true} [#sigma_{x}="+l+
             "];#Delta#alpha (mrad);Events").c_str(),
            200, -5, 5);
        // Vertex z residual vs true opening angle — shows where resolution fails
        H["ResZ_vsAlpha"   +t] = new TH2F(("ResZ_vsAlpha"   +t).c_str(),
            ("#Deltaz vs #alpha_{true} [#sigma_{x}="+l+
             "];#alpha_{true} (mrad);#Deltaz (mm)").c_str(),
            100, 0, 100, 400, -4000, 4000);
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
            if (particle.getPDG() == Cuts::LAMBDA_PDG &&
                particle.getParents().size() == 2)
            {
                double ez = particle.getEndpoint().z;
                h1("LambdaEnd")->Fill(ez);
                if      (ez <  6000)                h1("LambdaEnd_1")->Fill(ez);
                else if (ez >= 6000 && ez <= 20800) h1("LambdaEnd_2")->Fill(ez);
                else                                h1("LambdaEnd_3")->Fill(ez);

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
            // getDaughters().size()==0: Lambda escaped the detector boundary
            // or its decay was not simulated. Fill endpoint histogram and skip.
            if (particle.getDaughters().size() == 0) {
                h1("LambdaZnoDaughters")->Fill(particle.getEndpoint().z);
                continue;
            }

            // ── Geometric decay vertex cut ─────────────────────────────────
            // Accept only Lambdas decaying in the field-free region between the
            // last magnet (DECAY_Z_MIN = 20.8 m) and the ZDC face
            // (DECAY_Z_MAX = 35 m). In this region daughters travel in straight
            // lines, so the ZDC hit position is computed by simple geometry.
            // The daughter count cut (==2) selects Lambda -> p + pi-,
            // rejecting Lambda -> n + pi0 (also 2 daughters but neutral).
            double ez = particle.getEndpoint().z;
            if (ez <= Cuts::DECAY_Z_MIN ||
                ez >= Cuts::DECAY_Z_MAX ||
                (int)particle.getDaughters().size() != Cuts::LAMBDA_N_DAUGHTERS)
                continue;

            counters.lambda++;

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
                else if (pdg == Cuts::PDG_NEUTRON) {
                    DaughterKin kin_n;
                    kin_n.compute(mom.x, mom.y, mom.z, xv, yv, zv, Cuts::ZDC_Z, Cuts::CROSSING_ANGLE);
                    h1("Vertex_n")  ->Fill(daughter.getVertex().z);
                    h1("Endpoint_n")->Fill(daughter.getEndpoint().z);
                    h1("Theta_n")   ->Fill(kin_n.theta);
                    h2("XYatZDC_n") ->Fill(kin_n.xAtZDC, kin_n.yAtZDC);
                    counters.neutron++;
                }
                else if (pdg == Cuts::PDG_PION_ZERO) {
                    DaughterKin kin_pi0;
                    kin_pi0.compute(mom.x, mom.y, mom.z, xv, yv, zv, Cuts::ZDC_Z, Cuts::CROSSING_ANGLE);
                    h1("Vertex_pi0")  ->Fill(daughter.getVertex().z);
                    h1("Endpoint_pi0")->Fill(daughter.getEndpoint().z);
                    h1("Theta_pi0")   ->Fill(kin_pi0.theta);
                    h2("XYatZDC_pi0") ->Fill(kin_pi0.xAtZDC, kin_pi0.yAtZDC);
                    counters.pion_zero++;
                }
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

                // ── Analytical tracker resolution study ───────────────────
                // D = distance from decay point to first tracker plane [mm]
                //   |decay|...|plane1|---LENGTH---|plane4|---STANDOFF---|ZDC|
                double alpha_rad = alpha * M_PI / 180.0;  // rad
                double D = dist2ZDC - Tracker::STANDOFF - Tracker::LENGTH;

                // Guard: skip if decay is past plane1, or alpha too small
                if (D > 0.0 && alpha_rad > 1e-6) {
                    for (int i = 0; i < Tracker::N_SC; ++i) {
                        string t = Form("_sc%d", i);

                        double s_th    = sigma_theta_rad(Tracker::SIGMA_X[i]);
                        double s_alpha = std::sqrt(2.0) * s_th;         // rad
                        double s_z     = std::sqrt(2.0) * D * s_th / alpha_rad; // mm

                        h1("SigmaZ"         +t)->Fill(s_z);
                        h1("SigmaAlpha"     +t)->Fill(s_alpha * 1e3);       // mrad
                        h2("SigmaZ_vsAlpha" +t)->Fill(alpha_rad * 1e3, s_z);
                        h2("SigmaZ_vsD"     +t)->Fill(D, s_z);
                    }

                    // ── Smearing study ────────────────────────────────────
                    // Smear the true track angles with Gaussian noise whose
                    // width equals sigma_theta for each scenario, then
                    // reconstruct the opening angle and decay vertex z from
                    // the smeared directions.
                    //
                    // Track directions are represented as unit vectors.
                    // The tracker measures hits in x and y; fitting a straight
                    // line gives the direction (theta, phi). We smear theta
                    // and phi independently:
                    //   - theta smearing: sigma_theta (rad)
                    //   - phi smearing:   sigma_theta / sin(theta)
                    //     (a transverse displacement sigma_x maps to a larger
                    //      phi shift at small theta)
                    //
                    // Vertex reconstruction: each smeared track is a ray
                    // starting at plane1 (z = ZDC_Z - STANDOFF - LENGTH)
                    // in direction (theta_reco, phi_reco). The reconstructed
                    // vertex is the point of closest approach (PCA) of the
                    // two rays, extrapolated backward.
                    //
                    // PCA in the small-angle approximation:
                    //   The two rays at plane1 have positions (x1,y1) and
                    //   (x2,y2) and directions (tx1,ty1) and (tx2,ty2) where
                    //   tx = tan(theta)*cos(phi), ty = tan(theta)*sin(phi).
                    //   In the forward (large pz) limit the z at which they
                    //   cross in x is:
                    //     z_vtx = z_plane1 - (x2-x1)/(tx2-tx1)   [x-projection]
                    //   and similarly in y. We average the two projections.

                    // True angles in radians
                    double th_p_true  = kin_p.theta  * M_PI / 180.0;
                    double ph_p_true  = kin_p.phi    * M_PI / 180.0;
                    double th_pi_true = kin_pi.theta * M_PI / 180.0;
                    double ph_pi_true = kin_pi.phi   * M_PI / 180.0;

                    // True hit positions at plane1
                    // plane1 is at z_p1 = ZDC_Z - STANDOFF - LENGTH
                    double z_p1 = Cuts::ZDC_Z - Tracker::STANDOFF - Tracker::LENGTH;
                    double dz_p1_p  = z_p1 - zv;  // from decay to plane1
                    double dz_p1_pi = z_p1 - zv;  // same for both (same vertex)

                    double x_p_p1  = xv + dz_p1_p  * std::tan(th_p_true)  * std::cos(ph_p_true);
                    double y_p_p1  = yv + dz_p1_p  * std::tan(th_p_true)  * std::sin(ph_p_true);
                    double x_pi_p1 = xv + dz_p1_pi * std::tan(th_pi_true) * std::cos(ph_pi_true);
                    double y_pi_p1 = yv + dz_p1_pi * std::tan(th_pi_true) * std::sin(ph_pi_true);

                    for (int i = 0; i < Tracker::N_SC; ++i) {
                        string t = Form("_sc%d", i);

                        double s_th = sigma_theta_rad(Tracker::SIGMA_X[i]);

                        // Smear theta for each track
                        double th_p_reco  = th_p_true  + gRandom->Gaus(0.0, s_th);
                        double th_pi_reco = th_pi_true + gRandom->Gaus(0.0, s_th);

                        // Smear phi: sigma_phi = sigma_theta / sin(theta)
                        // Protect against sin(theta) ~ 0 (extremely forward tracks)
                        double s_phi_p  = (std::sin(th_p_true)  > 1e-6) ?
                                          s_th / std::sin(th_p_true)  : s_th;
                        double s_phi_pi = (std::sin(th_pi_true) > 1e-6) ?
                                          s_th / std::sin(th_pi_true) : s_th;

                        double ph_p_reco  = ph_p_true  + gRandom->Gaus(0.0, s_phi_p);
                        double ph_pi_reco = ph_pi_true + gRandom->Gaus(0.0, s_phi_pi);

                        // Reconstructed opening angle from smeared unit vectors
                        double dot_reco =
                            std::sin(th_p_reco)  * std::cos(ph_p_reco)  *
                            std::sin(th_pi_reco) * std::cos(ph_pi_reco)
                          + std::sin(th_p_reco)  * std::sin(ph_p_reco)  *
                            std::sin(th_pi_reco) * std::sin(ph_pi_reco)
                          + std::cos(th_p_reco)  * std::cos(th_pi_reco);
                        dot_reco = std::max(-1.0, std::min(1.0, dot_reco));
                        double alpha_reco = std::acos(dot_reco); // rad

                        // Reconstructed track slopes at plane1
                        // tx = tan(theta)*cos(phi), ty = tan(theta)*sin(phi)
                        double tx_p  = std::tan(th_p_reco)  * std::cos(ph_p_reco);
                        double ty_p  = std::tan(th_p_reco)  * std::sin(ph_p_reco);
                        double tx_pi = std::tan(th_pi_reco) * std::cos(ph_pi_reco);
                        double ty_pi = std::tan(th_pi_reco) * std::sin(ph_pi_reco);

                        // PCA vertex reconstruction: z where tracks cross in x and y
                        // Protect against parallel tracks (dtx or dty ~ 0)
                        double z_reco = zv; // fallback: true value
                        double dtx = tx_p - tx_pi;
                        double dty = ty_p - ty_pi;

                        if (std::abs(dtx) > 1e-9 && std::abs(dty) > 1e-9) {
                            double z_from_x = z_p1 - (x_p_p1 - x_pi_p1) / dtx;
                            double z_from_y = z_p1 - (y_p_p1 - y_pi_p1) / dty;
                            z_reco = 0.5 * (z_from_x + z_from_y); // average
                        } else if (std::abs(dtx) > 1e-9) {
                            z_reco = z_p1 - (x_p_p1 - x_pi_p1) / dtx;
                        } else if (std::abs(dty) > 1e-9) {
                            z_reco = z_p1 - (y_p_p1 - y_pi_p1) / dty;
                        }

                        // Residuals
                        double res_z     = z_reco - zv;               // mm
                        double res_alpha = (alpha_reco - alpha_rad) * 1e3; // mrad

                        h1("ResZ"         +t)->Fill(res_z);
                        h1("ResAlpha"     +t)->Fill(res_alpha);
                        h2("ResZ_vsAlpha" +t)->Fill(alpha_rad * 1e3, res_z);
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
        TCanvas* c = make_canvas("cLambdaEnd", "Lambda decay endpoint", 2, 2);
        c->cd(1); h1("LambdaEnd")  ->Draw();
        c->cd(2); h1("LambdaEnd_1")->Draw();
        c->cd(3); h1("LambdaEnd_2")->Draw();
        c->cd(4); h1("LambdaEnd_3")->Draw();
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
// ============================================================================
vector<TCanvas*> plot_tracker_study() {
    vector<TCanvas*> pages;

    const int colors[Tracker::N_SC] = {
        kBlue, kRed, kGreen+2, kOrange+7, kViolet
    };

    auto make_canvas = [&](const char* name, const char* title,
                           int nx = 1, int ny = 1) -> TCanvas* {
        TCanvas* c = new TCanvas(name, title, 800, 600);
        if (nx > 1 || ny > 1) c->Divide(nx, ny);
        pages.push_back(c);
        return c;
    };

    auto save_legend = [&](TLegend* leg) {
        leg->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i)
            leg->AddEntry(h1("SigmaZ" + string(Form("_sc%d",i))),
                          Tracker::LABELS[i], "l");
        leg->Draw();
    };

    // ── Analytical: page 1 — sigma_z overlay ─────────────────────────────
    {
        TCanvas* c = make_canvas("cSigmaZ", "Vertex z-resolution (analytical)");
        c->SetLogy();
        TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string t = Form("_sc%d", i);
            TH1F* h  = h1("SigmaZ" + t);
            h->SetLineColor(colors[i]);
            h->SetLineWidth(2);
            h->Draw(i == 0 ? "HIST" : "HIST SAME");
            leg->AddEntry(h, Tracker::LABELS[i], "l");
        }
        leg->Draw();
        c->Update();
    }

    // ── Analytical: page 2 — sigma_alpha overlay ─────────────────────────
    {
        TCanvas* c = make_canvas("cSigmaAlpha", "Opening-angle resolution (analytical)");
        c->SetLogy();
        TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string t = Form("_sc%d", i);
            TH1F* h  = h1("SigmaAlpha" + t);
            h->SetLineColor(colors[i]);
            h->SetLineWidth(2);
            h->Draw(i == 0 ? "HIST" : "HIST SAME");
            leg->AddEntry(h, Tracker::LABELS[i], "l");
        }
        leg->Draw();
        c->Update();
    }

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
    // For each scenario compute efficiency = fraction of events with
    // sigma_z < threshold, for each threshold in Tracker::EFF_THRESH.
    // Also compute mean sigma_z for the summary graph.
    {
        TCanvas* c = make_canvas("cSummary",
            "Mean #sigma_{z} and efficiency vs tracker pitch", 1, 2);

        std::vector<double> sx(Tracker::N_SC);
        std::vector<double> sy(Tracker::N_SC);
        std::vector<double> ey(Tracker::N_SC);

        // Efficiency storage: [scenario][threshold]
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
                    // Find the bin corresponding to the threshold
                    int bin = hSZ->FindBin(Tracker::EFF_THRESH[j]);
                    double below = hSZ->Integral(1, bin);
                    eff[i][j] = below / total;
                }
            }
        }

        // Top pad: mean sigma_z vs sigma_x
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

        // Bottom pad: efficiency curves, one per threshold
        c->cd(2);
        gPad->SetLogx();
        const int thr_colors[3] = {kBlue, kRed, kGreen+2};
        TLegend* legE = new TLegend(0.15, 0.15, 0.55, 0.45);
        legE->SetHeader("#sigma_{z} threshold", "C");
        for (int j = 0; j < Tracker::N_THRESH; ++j) {
            std::vector<double> effv(Tracker::N_SC);
            for (int i = 0; i < Tracker::N_SC; ++i) effv[i] = eff[i][j];
            TGraph* gEff = new TGraph(Tracker::N_SC, sx.data(), effv.data());
            gEff->SetTitle(
                "Efficiency vs tracker pitch"
                ";#sigma_{x} (mm);Fraction of events");
            gEff->SetMarkerStyle(20 + j);
            gEff->SetMarkerColor(thr_colors[j]);
            gEff->SetLineColor(thr_colors[j]);
            gEff->SetLineWidth(2);
            gEff->Draw(j == 0 ? "APL" : "PL SAME");
            gEff->GetYaxis()->SetRangeUser(0.0, 1.05);
            legE->AddEntry(gEff,
                Form("#sigma_{z} < %.0f mm", Tracker::EFF_THRESH[j]), "lp");
        }
        legE->Draw();
        c->Update();

        // Print summary table
        fmt::print("\n========== Tracker Resolution Study ==========\n");
        fmt::print("  {} planes, L={:.0f} mm, standoff={:.0f} mm from ZDC\n",
                   Tracker::N_PLANES, Tracker::LENGTH, Tracker::STANDOFF);
        fmt::print("  {:>12s}  {:>16s}  {:>16s}  {:>18s}",
                   "sigma_x (mm)", "sigma_th (urad)", "sigma_al (mrad)",
                   "mean sigma_z (mm)");
        for (int j = 0; j < Tracker::N_THRESH; ++j)
            fmt::print("  {:>14s}", Form("eff<%.0fmm", Tracker::EFF_THRESH[j]));
        fmt::print("\n");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            double st = sigma_theta_rad(Tracker::SIGMA_X[i]) * 1e6;
            double sa = std::sqrt(2.0) * sigma_theta_rad(Tracker::SIGMA_X[i]) * 1e3;
            fmt::print("  {:>12.3f}  {:>16.1f}  {:>16.4f}  {:>14.1f}+/-{:.1f}",
                       Tracker::SIGMA_X[i], st, sa, sy[i], ey[i]);
            for (int j = 0; j < Tracker::N_THRESH; ++j)
                fmt::print("  {:>14.1f}%", eff[i][j] * 100.0);
            fmt::print("\n");
        }
        fmt::print("==============================================\n\n");
    }

    // ── Smearing: page 6 — vertex z residuals overlay ────────────────────
    {
        TCanvas* c = make_canvas("cResZ",
            "Vertex z residual (smearing)", 1, 1);
        c->SetLogy();
        TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string t = Form("_sc%d", i);
            TH1F*  h = h1("ResZ" + t);
            h->SetLineColor(colors[i]);
            h->SetLineWidth(2);
            h->Draw(i == 0 ? "HIST" : "HIST SAME");
            leg->AddEntry(h, Tracker::LABELS[i], "l");
        }
        leg->Draw();
        c->Update();
    }

    // ── Smearing: page 7 — opening angle residuals overlay ───────────────
    {
        TCanvas* c = make_canvas("cResAlpha",
            "Opening-angle residual (smearing)", 1, 1);
        c->SetLogy();
        TLegend* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
        leg->SetHeader("Tracker pitch #sigma_{x}", "C");
        for (int i = 0; i < Tracker::N_SC; ++i) {
            string t = Form("_sc%d", i);
            TH1F*  h = h1("ResAlpha" + t);
            h->SetLineColor(colors[i]);
            h->SetLineWidth(2);
            h->Draw(i == 0 ? "HIST" : "HIST SAME");
            leg->AddEntry(h, Tracker::LABELS[i], "l");
        }
        leg->Draw();
        c->Update();
    }

    // ── Smearing: page 8 — residual vs opening angle, per scenario ───────
    {
        TCanvas* c = make_canvas("cResZ_vsAlpha",
            "#Deltaz vs opening angle (smearing)", 3, 2);
        for (int i = 0; i < Tracker::N_SC; ++i) {
            c->cd(i+1);
            h2("ResZ_vsAlpha" + string(Form("_sc%d",i)))->Draw("COLZ");
        }
        c->Update();
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

    vector<TCanvas*> pages = plot();
    save_plots(pages, Cuts::OUTPUT_PDF);

    vector<TCanvas*> tracker_pages = plot_tracker_study();
    save_tracker_study(tracker_pages);
}



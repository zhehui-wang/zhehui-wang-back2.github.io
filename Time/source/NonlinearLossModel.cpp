#include "NonlinearLossModel.h"

/*
 * References
 * [1] M. Dinu, F. Quochi, and H. Garcia, “Third-order nonlinearities in silicon at telecom wavelengths,” Appl. Phys. Lett., vol.
 * 82, no. 18, pp. 2954–2956, Apr. 2003. [2] K. Dietmar, Silicon-Organic Hybrid Platform for Photonic Integrated Circuits. KIT
 * Scientific Publishing, 2015.
 */

class NonlinearLossModel;

#ifdef __DEBUG_
NonlinearLossModel::NonlinearLossModel()
    : kPhotonEnergy1550(1.28e-19)
    , kStartPosition(0)
    , kStepLength(1e-7)
    , kWavelength(1550) {
    // default constructor should only be used for bebugging.
    initParameters();
}
#endif

double utils::loss_dB2scalar(double loss_in_dB) { return pow(10, -loss_in_dB / 10); }
double utils::loss_scaler2dB(double loss_in_scalar) { return -10 * log10(loss_in_scalar); }

NonlinearLossModel::NonlinearLossModel(double coupled_in_power, double waveguideLength)
    : kPhotonEnergy1550(1.28e-19)
    , kStartPosition(0)
    , kStepLength(1e-7)
    , kWavelength(laser_wavelength) {

    setParametersFromFile();
    setInitialPower(coupled_in_power);
    setWaveguideLength(waveguideLength);
}

void NonlinearLossModel::initParameters() {
    alpha_dB_ = 0.2;             // Default: Linear loss: 0.2dB/cm
    carrier_lifetime_ = 4e-9;    // Default: 4e-9s
    effective_mode_area_ = 1e-8; // Default: 1e-8cm^2 (1um^2)
    stop_position_ = 1.5;        // Length of the wavelength in question. Default: 1.5cm
    initial_power_ = 100e-3;     // The coupled power at distance 0. Default: 100e-3W (100mW)

    beta_ = 8e-10; // Default: TPA coefficient at 1550nm [1]
    sigma_ = 1.45e-17 * pow(kWavelength / 1550,
                            2); // Default: FCA coefficient: 1.45*10^17 cm2 [2]. sigma_ = 1.45e-17*pow(kWavelength / 1550, 2);

    updateModel();
}

// Invoke updateModel() whenever a parameter is changed
void NonlinearLossModel::updateModel() {
    // model parameters generation
    photon_energy_ = kPhotonEnergy1550 * (1550 / kWavelength); // Default: hv at 1.55um is 1.28 * 10^-19 J
    alpha_ = alpha_dB_ * log(10) / 10;
    gamma_ = (sigma_ * carrier_lifetime_ * beta_) / (2 * photon_energy_);
    initial_intensity_ = initial_power_ / effective_mode_area_;

    // reset the buffered results
    totalLoss_ = -1;
}

void NonlinearLossModel::setLinearLoss(double linearLoss_dB) {
    alpha_dB_ = linearLoss_dB;
    updateModel();
}
void NonlinearLossModel::setCarrierLifetime(double lifetime) {
    carrier_lifetime_ = lifetime;
    updateModel();
}
void NonlinearLossModel::setEffectiveModeArea(double Aeff) {
    effective_mode_area_ = Aeff;
    updateModel();
}
void NonlinearLossModel::setWaveguideLength(double length) {
    stop_position_ = length;
    updateModel();
}
void NonlinearLossModel::setInitialPower(double power) {
    initial_power_ = power;
    updateModel();
}

// loaders
void NonlinearLossModel::resetParameters() {
    initParameters();
    updateModel();
}

void NonlinearLossModel::setParametersFromFile() {
    alpha_dB_ = propagation_loss * 10 / log(10); // Default: Linear loss: 0.2dB/cm
    carrier_lifetime_ = carrier_lifetime;        // Default: 4e-9s
    effective_mode_area_ = effective_mode_area;  // Default: 1e-8cm^2 (1um^2)

    beta_ = TPA_coefficient;  // Default: TPA coefficient at 1550nm [1]
    sigma_ = FCA_coefficient; // Default: FCA coefficient: 1.45*10^17 cm2 [2]. sigma_ = 1.45e-17 * pow(wavelength / 1.55, 2);
    updateModel();
}

void NonlinearLossModel::echoParameters() {
    using namespace std;
    cout << "------------------Start of Nonlinear Model Parameters-----------------------" << endl;
    //  (Only C-band is supported so far)
    cout << "Central wavelength: " << kWavelength << "nm" << endl;
    cout << "Linear loss per centimeter: " << alpha_dB_ << "dB/cm" << endl;
    cout << "Effective mode area: " << effective_mode_area_ * 1e8 << "um^2" << endl;
    cout << "Carrier lifetime : " << carrier_lifetime_ * 1e9 << "ns" << endl;
    cout << "Waveguide length: " << stop_position_ << "cm" << endl;
    cout << "Two-photon absorption (TPA) coefficient: " << beta_ << "cm/W" << endl;
    cout << "Free-carrier absorption (FCA) coefficient: " << sigma_ << "cm^2" << endl;
    cout << "Coupled-in power: " << initial_power_ << "W" << endl;
    cout << "Initla intensity: " << initial_intensity_ << "W/cm^2" << endl;
    cout << "--------------------End of Nonlinear Model Parameters-----------------------" << endl;
}

void NonlinearLossModel::operator()(const intensity_state& I, intensity_state& dIdz, const double z) const {
    // dIdz = -alpha_ * I - beta_ * pow(I, 2) - gamma_ * pow(I, 3);
    dIdz[0] = -alpha_ * I[0] - beta_ * pow(I[0], 2) - gamma_ * pow(I[0], 3);
}

double NonlinearLossModel::totalLossAcrossWaveguide() {
    if (totalLoss_ != -1) {
        return totalLoss_;
    }

    // Calculate the loss if no buffered result is found
    intensity_state intensity{initial_intensity_};
    boost::numeric::odeint::bulirsch_stoer<intensity_state> stepper;

    using namespace std;
    cout << "start_position: " << kStartPosition << "cm stop_position_: " << stop_position_ << "cm" << endl;
    cout << "initial power: " << initial_power_ << "W initial intensity: " << intensity[0] << "W/cm^2" << endl;
    cout << "processing....................."
         << "\n.............."
         << "\n........" << endl;
    size_t steps = boost::numeric::odeint::integrate_const(stepper, *this, intensity, kStartPosition, stop_position_, 1e-7);

    totalLoss_ = 10 * log10(initial_intensity_ / intensity[0]);
    cout << "final power: " << intensity[0] * effective_mode_area_ << "W final intensity: " << intensity[0] << "W/cm^2" << endl;
    cout << "total loss(dB): " << totalLoss_ << endl;
    return totalLoss_;
}

/**
 */
double NonlinearLossModel::totalLossAcrossWaveguide(double initalPower) {
    setInitialPower(initalPower);
    return totalLossAcrossWaveguide();
}

double NonlinearLossModel::nonlinearLossAcrossWaveguide() { return totalLossAcrossWaveguide() - linearLossAcrossWaveguide(); }

double NonlinearLossModel::nonlinearLossAcrossWaveguide(double initalPower) {
    return totalLossAcrossWaveguide(initalPower) - linearLossAcrossWaveguide();
}

double NonlinearLossModel::linearLossAcrossWaveguide() { return (stop_position_ - kStartPosition) * alpha_dB_; }

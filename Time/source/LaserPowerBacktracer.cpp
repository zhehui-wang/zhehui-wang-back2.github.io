#include "LaserPowerBacktracer.h"

class LaserPowerBacktracer;

#ifdef __DEBUG_
/**
 * @brief The default construct of LaserPowerBacktracer class should only be used for debugging purpose.
 *
 */
LaserPowerBacktracer::LaserPowerBacktracer()
    : NonlinearLossModel() {}
#endif

/**
 * @brief Construct a new LaserPowerBacktracer object by specifying the coupled-out power (in Watt) by the end of the waveguide,
 * and the waveguide length.
 *
 * @param powerInWattBeforeReceiver The coupled-out power by the end of the waveguide.
 * @param waveguideLength Length of the waveguide in question.
 */
LaserPowerBacktracer::LaserPowerBacktracer(double powerInWattBeforeReceiver, double waveguideLength)
    : NonlinearLossModel(-1, waveguideLength)
    , final_power_(powerInWattBeforeReceiver) {
    updateModel();
}

void LaserPowerBacktracer::echoParameters() const {
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
    // cout << "Coupled-out power: " << final_power_ * 1e3 << "mW" << endl;
    // cout << "Final intensity: " << final_intensity_ << "W/cm^2" << endl;
    cout << "--------------------End of Nonlinear Model Parameters-----------------------" << endl;
}

double LaserPowerBacktracer::totalLossAcrossWaveguide() {
    if (totalLoss_ != -1) {
        return totalLoss_;
    }

    // Calculate the loss if no buffered result is found
    intensity_state intensity{final_intensity_};
    boost::numeric::odeint::bulirsch_stoer<intensity_state> stepper;

    using namespace std;

#ifdef __DEBUG_
    cout << "\nstart_position: " << getStartPosition() << " stop_position: " << getStopPosition() << "(cm)" << endl;
    cout << "final power: " << final_power_ << "(W) final intensity: " << final_intensity_ << "(W/cm^2)" << endl;
    cout << "processing....................."
         << "\n.............."
         << "\n........" << endl;
#endif
    size_t steps =
        boost::numeric::odeint::integrate_const(stepper, *this, intensity, getStartPosition(), getStopPosition(), 1e-7);

    totalLoss_ = 10 * log10(intensity[0] / final_intensity_);
#ifdef __DEBUG_
    cout << "initial power: " << intensity[0] * effective_mode_area_ << " initial intensity: " << intensity[0] << endl;
    cout << "total loss(dB): \n" << totalLoss_ << endl;
#endif
    return totalLoss_;
}

double LaserPowerBacktracer::nonlinearLossAcrossWaveguide() {
    return LaserPowerBacktracer::totalLossAcrossWaveguide() - linearLossAcrossWaveguide();
}

/**
 * @brief This function sets the required power level by the end of the waveguide. The nonlinear loss will the be calculated by
 * backtracing the coupled-in power.
 *
 * @param finalPower The required power level in Watt by the end of the waveguide.
 */
void LaserPowerBacktracer::setFinalPower(double finalPower) {
    final_power_ = finalPower;
    updateModel();
}

/**
 * @brief This function intends to recalculate the model-related parameters from the user-input parameters. This function should
 * be invoked at leaser once before nonlinear loss is calculated.
 */
void LaserPowerBacktracer::updateModel() {
    final_intensity_ = final_power_ / effective_mode_area_;
    totalLoss_ = -1; // Reset buffer
}

/**
 * @brief This is a callback function for integating the optical power intensity over the wavelength, z. See
 * boost::numeric::odeint for details.
 *
 * @param I denotes the optical power intensity at a given position on the waveguide. (Unit: Watt)
 * @param dIdz denotes the changing rate of power intensity at a given position.
 * @param z denote the postion on the given waveguide.
 */
void LaserPowerBacktracer::operator()(const intensity_state& I, intensity_state& dIdz, const double z) const {
    // This is the reverse integral of the nonlinear loss model, which derive final power from the initial coupled-in power
    // In other words, this function backtraces the initial power given the final power at the end of the waveguide
    // dIdz = -alpha * I - beta * pow(I, 2) - gamma * pow(I, 3);
    dIdz[0] = +alpha_ * I[0] + beta_ * pow(I[0], 2) + gamma_ * pow(I[0], 3);
}
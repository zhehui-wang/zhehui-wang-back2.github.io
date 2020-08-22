/*********************************************************************************
 *
 * File name:			OEIL.cpp
 * Version:			3.0
 * Software:
 * Authors:			Zhehui Wang, Jiang Xu
 * Website:			http://www.ece.ust.hk/~eexu/
 * The copyright information of this program can be found in the file COPYRIGHT.
 *
 *********************************************************************************/
#include "definition_electrical.h"
#include "definition_optical.h"

#include "LaserPowerBacktracer.h"
#include <stdio.h>

/*the following value is for the constants*/
#define PI 3.1415926
#define E 2.7182818

/*the following value is for the values that are estimated coefficient*/
#define DELTAV 0.1     // the saturation voltage of the transisitor unit:V
#define LIGHTSPEED 30  // the light speed in vacuum  unitcm/ns
#define ILADENSITY 0.3 // the working current of the limiting amplifier for each GHz unit:mA

/*sensitivity is used to calculate the optical_modulation_amplitude of receiver*/
void sensitivity_optical() {
    sensitivity_oma = (0.0001 * tia_noise_density * signal_to_noise_ratio * pow(data_rate_optical / 20, 0.5) +
                       la_voltage_threshold * 0.002 / tia_transimpendance) /
                      pd_responsity; // the total sensitivity of optical modulation amplitude
}

/*micorresonator is used to calculate the real radius of microresonator, the wavelength spacing, and the split factor r*/
void microresonator() {
    double fsr, phase;
    phase = 4 * PI * PI * mr_refractive_index * mr_radius_range / laser_wavelength * 1000; // the estimated phase
    phase = (int)(phase / (2 * PI)) * 2 * PI; // the real phase at resonace wavelength
    mr_radius_real = phase * laser_wavelength / (1000 * 4 * PI * PI * mr_refractive_index); // the real radius of the mr
    fsr = laser_wavelength * laser_wavelength /
          (2 * PI * mr_refractive_index * mr_radius_real * 1000); // free spectrum range of MR, the unit: nm
    delta_wavelength = fsr / (double)number_of_wavelengths;       // the spacing betweeen two naerby wavelength
    mr_power_split_r = sqrt(1 - pow(mr_power_split_k, 2.0));      // the relationship between r and k is that r^2+k^2=1
}

/*crosstalk is used to calculate the crosstalk noise coefficiencts in WDM optical interconnects*/
void crosstalk_optical() {
    double phase, drop_transmission_factor;
    for (int i = 1; i <= number_of_wavelengths / 2; i++) {
        phase = 4 * PI * PI * mr_refractive_index * mr_radius_real / (laser_wavelength + i * delta_wavelength) *
                1000; // the phase of the signal with given wavelength
        drop_transmission_factor = pow((1 - pow(mr_power_split_r, 2.0)), 2.0) *
                                   mr_attenuation // the transmission ratio on the drop port
                                   / (1 - 2 * pow(mr_power_split_r, 2.0) * mr_attenuation * cos(phase) +
                                      pow(mr_power_split_r, 4.0) * pow(mr_attenuation, 2.0));
        crosstalk_coefficient_optical += 2 * drop_transmission_factor; // the total crosstalk noise coefficient
    }
}

/*attenuation is used to calculate the total attenuation of the entire optical interconnect*/
void attenuation_optical() {
    double phase, pass_transmission_factor, loss_insertion;
    total_attenuation_optical = 1;
    for (int i = 1; i < number_of_wavelengths; i++) {
        phase = 4 * PI * PI * mr_refractive_index * mr_radius_real / (laser_wavelength - i * delta_wavelength) *
                1000; // the phase of the signal with given wavelength
        pass_transmission_factor =
            (pow(mr_power_split_r, 2.0) * pow(mr_attenuation, 2.0) -
             2 * pow(mr_power_split_r, 2.0) * mr_attenuation * cos(phase) + pow(mr_power_split_r, 2.0)) /
            (1 - 2 * pow(mr_power_split_r, 2.0) * mr_attenuation * cos(phase) +
             pow(mr_power_split_r, 4.0) * pow(mr_attenuation, 2.0)); // the transmission ratio on the pass port
        total_attenuation_optical *= pass_transmission_factor;       // the total passing loss of mr
    }
    // phase = 4 * PI * PI * mr_refractive_index * mr_radius_real / (laser_wavelength)*1000; // the phase at resonace wavelength
    loss_insertion =
        pow((1 - pow(mr_power_split_r, 2.0)), 2.0) * mr_attenuation // the insertion loss of mr
        / (1 - 2 * pow(mr_power_split_r, 2.0) * mr_attenuation * 1 + pow(mr_power_split_r, 4.0) * pow(mr_attenuation, 2.0));

    if (!is_nonlinear_model_enabled) {
        total_attenuation_optical = pow(total_attenuation_optical, 2) * pow(optical_pin_loss, 2.0) *
                                    pow(E, -propagation_loss * length_optical) *
                                    loss_insertion; // calculate the total attenuation
        return;
    }
    double loss_passing = total_attenuation_optical;

    double P0, P1, P2, P3; // Intermediate powers
    double A1, A2, A3;     // Intermediate attenuations
    double total_receive_in_watt = sensitivity_oma * number_of_wavelengths * 1e-3;
    sensitivity_optical();
    P3 = total_receive_in_watt / (loss_insertion * loss_passing);

    LaserPowerBacktracer backtracer3(P3, coupler2receiver_distance);
    // Example of viewing the nonlinear model under the hood
    backtracer3.echoParameters();
    // std::cout << "The nonlinear loss across this segment is: " << backtracer3.nonlinearLossAcrossWaveguide() << "dB" <<
    // std::endl;
    A3 = utils::loss_dB2scalar(backtracer3.totalLossAcrossWaveguide());
    P2 = P3 / (A3 * pow(optical_pin_loss, 2.0) * pow(E, -propagation_loss * length_optical));

    LaserPowerBacktracer backtracer2(P2, modular2coupler_distance);
    A2 = utils::loss_dB2scalar(backtracer2.totalLossAcrossWaveguide());
    P1 = P2 / (A2 * loss_passing);

    LaserPowerBacktracer backtracer1(P1, laser2modular_distance);
    A1 = utils::loss_dB2scalar(backtracer1.totalLossAcrossWaveguide());
    P0 = P1 / (A1 * loss_insertion * loss_passing);
    // printf("P0: %f mW\n", P0 * 1e3);

    std::cout << "The total nonlinear loss is: "
              << backtracer1.totalLossAcrossWaveguide() + backtracer2.totalLossAcrossWaveguide() +
                     backtracer3.totalLossAcrossWaveguide()
              << "dB\n"
              << std::endl;
    total_attenuation_optical = A1 * A2 * A3 * pow(loss_passing, 2) * pow(optical_pin_loss, 2.0) * loss_insertion *
                                pow(E, -propagation_loss * length_optical);
}

/*energy is used to calculate the energy consumption of the optical interconenct*/
void energy_optical() {
    double energy_optical = 0;   // power consumption of optical interconnect
    double energy_modulator = 0; // power consumption of modulator
    double energy_serdes = 0;    // power consumption of serdes
    energy_optical = (sensitivity_oma / (total_attenuation_optical *
                                         (1 - crosstalk_coefficient_optical - laser_extinction_ratio) * laser_slope_efficiency) +
                      laser_threshold_current) *
                         laser_voltage +
                     (0.001 * PI * data_rate_optical * pd_capacitance * DELTAV / 2 + ILADENSITY * data_rate_optical / 2) *
                         driver_voltage; // the total energy of the optical interconnect
    energy_modulator = mr_tuning_power + mr_static_power / 2 + mr_dynamic_power * data_rate_optical / 4;
    energy_serdes =
        9 * ceil(log((double)serdes_ratio_optical) / log(2)) * serdes_cur_optical * driver_voltage * data_rate_optical;

    energy_consumption_optical = (energy_optical + (1 - is_direct_modulation) * energy_modulator + energy_serdes) /
                                 data_rate_optical; // the energy consumed per bit
}

/*bandwidth is used to calculate the area bandwidth density and linear bandwdith density*/
void bandwidth_optical() {
    area_density_optical = 1000000 * number_of_wavelengths * data_rate_optical /
                           (optical_pin_height * optical_pin_width); // the area bandwidth density of optical interconnect
    linear_density_optical =
        1000 * number_of_wavelengths * data_rate_optical / wg_pitch; // the linear bandwidth density of optical interconnect
}

/*area is used to calculate the area size in optical interconnect*/
void area_size_optical() {
    double area_circuit = 0;          // area of electrical circuits.
    double area_optical_modulate = 0; // area of lasers and modulators.
    area_circuit = 9 * ceil(log((double)serdes_ratio_optical) / log(2)) * serdes_area_optical * data_rate_optical;
    area_optical_modulate = mr_area * 2 + laser_area; // two MRs and one laser
    area_optical = (area_circuit + area_optical_modulate) / 1000000;
}

/*latency is used to calculate the latency of optical interconnect*/
void latency_optical() {
    latency_value_optical = wg_refractive_index * length_optical / LIGHTSPEED +
                            (2 * serdes_ratio_optical - 1) / data_rate_optical; // the latency of optcal signals in waveguide
}

/*sensitivity is used to calculate the minimum required voltage of receiver*/
void sensitivity_electrical() { sensitivity_la = la_threshold_voltage; }

/*crosstalk is used to calculate the crosstalk noise coefficiencts in PCB traces*/
void crosstalk_electrical() {
    crosstalk_coefficient_electrical = 0;
    for (int i = 1; i < (number_of_parallel_traces / 2) + 1; i++) {
        crosstalk_coefficient_electrical +=
            2 * pow(pcb_layer_height, 2.0) /
                (4 * pow((i * pcb_trace_pair_pitch - 2 * pcb_trace_width), 2.0) + pow(pcb_layer_height, 2.0)) -
            4 * pow(pcb_layer_height, 2.0) / (4 * pow(i * pcb_trace_pair_pitch, 2.0) + pow(pcb_layer_height, 2.0)) +
            2 * pow(pcb_layer_height, 2.0) /
                (4 * pow((i * pcb_trace_pair_pitch + 2 * pcb_trace_width), 2.0) + pow(pcb_layer_height, 2.0));
    }
}

/*attenuation is used to calculate the total attenuation of the entire electrical interconnect*/
void attenuation_electrical() {
    double alpha =
        ((pcb_trace_width + pcb_trace_height) * trace_direct_current_r * pow((data_rate_electrical * 500 / 13.8), 0.5)) /
            (2 * trace_characteristic_z * pcb_trace_width) +
        PI * data_rate_electrical * 0.5 * trace_direct_current_r * 0.001 * pcb_loss_tangent * trace_characteristic_z;
    double pin_loss = 1 - pow(E, (-1 / (trace_characteristic_z * data_rate_electrical * electrical_pin_load_c * 0.001)));
    total_attenuation_electrical = pow(E, (-alpha * length_electrical)) * pow((pin_loss), 2);
}

/*energy is used to calculate the energy consumption of the electrical interconnect*/
void energy_electrical() {
    double energy_electrical = 0; // power consumption of electrical interconnect.
    double energy_serdes = 0;     // power consumption of serdes
    energy_electrical =
        circuit_voltage *
        (2 * (2 * sensitivity_la) /
             ((total_attenuation_electrical - crosstalk_coefficient_electrical - la_offset_coefficent) * trace_input_impendance) +
         ILADENSITY * data_rate_electrical / 2);
    energy_serdes =
        9 * ceil(log((double)serdes_ratio_electrical) / log(2)) * serdes_cur_electrical * circuit_voltage * data_rate_electrical;
    energy_consumption_electrical = (energy_electrical + energy_serdes) / data_rate_electrical; // the energy consumed per bit
}

/*bandwidth is used to calculate the area bandwidth density and linear bandwidth density*/
void bandwidth_electrical() {
    double data_rate_step = 0.1;
    double data_rate = 0;
    double margin = 1;
    double alpha, pin_loss;
    // loop and find the maximum acceptable data rate
    do {
        data_rate += data_rate_step;
        alpha = ((pcb_trace_width + pcb_trace_height) * trace_direct_current_r * pow((data_rate * 500 / 13.8), 0.5)) /
                    (2 * trace_characteristic_z * pcb_trace_width) +
                PI * data_rate * 0.5 * trace_direct_current_r * 0.001 * pcb_loss_tangent * trace_characteristic_z;
        pin_loss = 1 - pow(E, (-1 / (trace_characteristic_z * data_rate * electrical_pin_load_c * 0.001)));
        margin =
            pow(E, (-alpha * length_electrical)) * pow((pin_loss), 2) - crosstalk_coefficient_electrical - la_offset_coefficent;
    } while (margin > 0.01);
    // calculating the densities
    area_density_electrical = data_rate / 2 * package_pin_pitch * package_pin_pitch;
    linear_density_electrical = data_rate / (pcb_trace_pair_pitch * 0.0254);
}

/*area is used to calculate the area size in electrical interconnect*/
void area_size_electrical() {
    double area_circuit = 0; // area of electrical circuits.area_electrical=2;
    area_circuit = 9 * ceil(log((double)serdes_ratio_electrical) / log(2)) * serdes_area_electrical * data_rate_electrical;
    area_electrical = (area_circuit) / 1000000;
}

/*latency is used to calculate the latency of electrical interconnect*/
void latency_electrical() {
    latency_value_electrical =
        sqrt(pcb_dielectric) * length_electrical / LIGHTSPEED +
        (2 * serdes_ratio_electrical - 1) / data_rate_electrical; // the latency of electrical signals in trace
}

/*this is to make the program maintain hold position*/
void program_wait() {
    printf("***************      Program Ended       **************\n");
    exit(1);
}

/*main body*/
int main() {
    /*print head information*/
    printf("\n");
    printf("*******************************************************");
    printf("\n* Software:        OEIL 3.0                           *");
    printf("\n* Authors:         Zhehui Wang, Jiang Xu              *");
    printf("\n* Website:         http://www.ece.ust.hk/~eexu        *");
    printf("\n*******************************************************");
    printf("\n");
    printf("\n");

    // while (true) {
    //     double testNum;
    //     char name[150];
    //     char unit[50];
    //     scanf("%lf %s %s", &testNum, name, unit);
    //     printf("Here I can print out %lf %s %s\n", testNum, name, unit);
    //     if (testNum > -1.5 && testNum < -0.5)
    //         break;
    // }
    // return 0;

    /* input files for optical &electrical*/
    printf("**************      Configurations      ***************\n");
    printf("\n");
    read_configuration_optical();
    read_configuration_electrical();

    printf("**************        Parameters        ***************\n");
    printf("\n");
    read_parameter_optical();
    read_parameter_electrical();

    /* calculating for optical*/
    sensitivity_optical();
    microresonator();
    crosstalk_optical();
    attenuation_optical();
    energy_optical();
    bandwidth_optical();
    area_size_optical();
    latency_optical();

    /* calculating for electrical*/
    sensitivity_electrical();
    crosstalk_electrical();
    attenuation_electrical();
    energy_electrical();
    bandwidth_electrical();
    area_size_electrical();
    latency_electrical();

    /* output file*/
    printf("**************    Calculated  Results    **************\n");
    printf("\n");
    write_result_optical();
    write_result_electrical();
    program_wait();
}

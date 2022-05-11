#include "CellICCBioPhy.hpp"
#include <cmath>
#include <cassert>
#include <memory>
#include "Exception.hpp"
#include "OdeSystemInformation.hpp"
#include "RegularStimulus.hpp"
#include "HeartConfig.hpp"
#include "IsNan.hpp"

#include "OutputFileHandler.hpp"

using namespace std;

    double CellICCBioPhy::Get_ICC_Membrane__Cm()
    {
        return var_ICC_Membrane__Cm;
    }

    double CellICCBioPhy::Get_chaste_interface__i_ionic()
    {
        return var_chaste_interface__i_ionic;
    }

    double CellICCBioPhy::GetIntracellularCalciumConcentration()
    {
        // std::cout<< "Ca: " << var_chaste_interface__ICC_Membrane__Ca_i << std::endl;
        // return var_chaste_interface__ICC_Membrane__Ca_i;
        return mStateVariables[1];
    }

    CellICCBioPhy::CellICCBioPhy(boost::shared_ptr<AbstractIvpOdeSolver> pSolver, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCardiacCell(
                pSolver,
                22,
                0,
                pIntracellularStimulus)
    {
        // Time units: millisecond
        //
        this->mpSystemInfo = OdeSystemInformation<CellICCBioPhy>::Instance();
        Init();

        this->mParameters[0] = 0.0006; // (c,PU_unit__IP3) [millimolar]


		/******************************/
		/* FSM parameters initialized */
		/******************************/
        this->mParameters[1] /*V_excitation*/ = -67;
        this->mParameters[2] /*is_passed_threshold*/ = 0 /*false*/;
        this->mParameters[3] = 0;
        this->mParameters[4] = 0;
        this->mParameters[5] /*inactive_duration*/ = 0;
        this->mParameters[6] /*non-refractory period*/ = 7000;
        this->mParameters[7] /*ode-time-step*/ = 0.01;
        this->mParameters[8] /*is_first_cycle_passed*/ = 0 /*false*/;
        this->mParameters[9] /*t_start*/ = 5000;
        this->mParameters[10] /*total_time_elapsed*/ = 0;
        this->mParameters[11] /*is_negative_stroke_passed*/ = 0;
		/******************************/
		/* END */
		/******************************/

    }

    CellICCBioPhy::~CellICCBioPhy()
    {
    }

    void CellICCBioPhy::VerifyStateVariables()
    {

    }

    double CellICCBioPhy::GetIIonic(const std::vector<double>* pStateVariables)
    {
        // For state variable interpolation (SVI) we read in interpolated state variables,
        // otherwise for ionic current interpolation (ICI) we use the state variables of this model (node).
        if (!pStateVariables) pStateVariables = &rGetStateVariables();
        const std::vector<double>& rY = *pStateVariables;
        double var_chaste_interface__ICC_Membrane__Vm = (mSetVoltageDerivativeToZero ? this->mFixedVoltage : rY[0]);


		/******************************/
		/* FSM modifications */
		/******************************/

        if (mParameters[2] /*is_passed_threshold*/ < 50)
	    {
            return 0.0;
        }
		/******************************/
		/* END */
		/******************************/

        // Units: millivolt; Initial value: -67
        double var_chaste_interface__ICC_Membrane__Ca_i = rY[1];
        // Units: millimolar; Initial value: 1.06669634258941e-5
        double var_chaste_interface__d_Ltype__d_Ltype = rY[2];
        // Units: dimensionless; Initial value: 9.57461476818094e-6
        double var_chaste_interface__f_Ltype__f_Ltype = rY[3];
        // Units: dimensionless; Initial value: 0.982524908182641
        double var_chaste_interface__f_ca_Ltype__f_ca_Ltype = rY[4];
        // Units: dimensionless; Initial value: 0.999999999912651
        double var_chaste_interface__d_VDDR__d_VDDR = rY[5];
        // Units: dimensionless; Initial value: 0.00111781215370465
        double var_chaste_interface__f_VDDR__f_VDDR = rY[6];
        // Units: dimensionless; Initial value: 0.769965220525522
        double var_chaste_interface__d_CaCl__d_CaCl = rY[7];
        // Units: dimensionless; Initial value: 0.000362322744959391
        double var_chaste_interface__d_ERG__d_ERG = rY[8];
        // Units: dimensionless; Initial value: 0.199999791554489
        double var_chaste_interface__d_kv11__d_kv11 = rY[9];
        // Units: dimensionless; Initial value: 0.00441922841884097
        double var_chaste_interface__f_kv11__f_kv11 = rY[10];
        // Units: dimensionless; Initial value: 0.996587027588045
        double var_chaste_interface__d_Na__d_Na = rY[11];
        // Units: dimensionless; Initial value: 0.0161886058733976
        double var_chaste_interface__f_Na__f_Na = rY[12];
        // Units: dimensionless; Initial value: 0.166158939491521
        double var_chaste_interface__d_NSCC__d_NSCC = rY[13];
        // Units: dimensionless; Initial value: 0.00167688179762435

        var_ICC_Membrane__Cm = 0.025; // capacitance_units
        const double var_I_Na__I_Na = 20.0 * var_chaste_interface__f_Na__f_Na * var_chaste_interface__d_Na__d_Na * (var_chaste_interface__ICC_Membrane__Vm - 40.5723805548); // current_units
        const double var_I_Ltype__I_Ltype = 2.0 * var_chaste_interface__f_Ltype__f_Ltype * var_chaste_interface__d_Ltype__d_Ltype * var_chaste_interface__f_ca_Ltype__f_ca_Ltype * (var_chaste_interface__ICC_Membrane__Vm - (log(2.5 / var_chaste_interface__ICC_Membrane__Ca_i) * 13.3568673135)); // current_units
        const double var_I_VDDR__I_VDDR = 3.0 * var_chaste_interface__f_VDDR__f_VDDR * var_chaste_interface__d_VDDR__d_VDDR * (var_chaste_interface__ICC_Membrane__Vm - (log(2.5 / var_chaste_interface__ICC_Membrane__Ca_i) * 13.3568673135)); // current_units
        const double var_I_kv11__I_kv11 = 6.3 * var_chaste_interface__f_kv11__f_kv11 * var_chaste_interface__d_kv11__d_kv11 * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
        const double var_I_BK__I_BK = 37.3 * (1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm *  -0.0588235294118) - (2.0 * log(var_chaste_interface__ICC_Membrane__Ca_i * 1000.0))))) * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
        const double var_I_ERG__I_ERG = 2.5 * var_chaste_interface__d_ERG__d_ERG * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
        const double var_I_CaCl__I_CaCl = 10.1 * var_chaste_interface__d_CaCl__d_CaCl * (var_chaste_interface__ICC_Membrane__Vm -  -11.2332051638); // current_units
        const double var_I_NSCC__I_NSCC = 12.15 * var_chaste_interface__d_NSCC__d_NSCC * (var_chaste_interface__ICC_Membrane__Vm - 4.40291009976e-06); // current_units
        const double var_I_bk__I_bk = 0.15 * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
        var_chaste_interface__i_ionic = 0.001 * (((var_I_Na__I_Na + var_I_Ltype__I_Ltype + var_I_VDDR__I_VDDR + var_I_kv11__I_kv11 + var_I_ERG__I_ERG + var_I_BK__I_BK + var_I_CaCl__I_CaCl + var_I_NSCC__I_NSCC + var_I_bk__I_bk) / var_ICC_Membrane__Cm) * HeartConfig::Instance()->GetCapacitance()); // uA_per_cm2

        const double i_ionic = var_chaste_interface__i_ionic;
        EXCEPT_IF_NOT(!std::isnan(i_ionic));
        return i_ionic;
    }

    void CellICCBioPhy::EvaluateYDerivatives(double var_chaste_interface__Time__time, const std::vector<double>& rY, std::vector<double>& rDY)
    {


		/******************************/
		/* FSM modifications */
		/******************************/

        double ode_time_step = mParameters[7] /*ode-time-step*/;
        mParameters[10] /*total_time_elapsed*/ += ode_time_step;
        if (( rY[0] > mParameters[1] /*V_excitation*/ ) || ( mParameters[5] /*inactive_duration*/ >= mParameters[6] /*non-refractory period*/ ) || ((mParameters[10] /*total_time_elapsed*/>=mParameters[9] /*t_start*/) && (mParameters[8] /*is_first_cycle_passed*/ < 100)))
        {
            mParameters[2] /*is_passed_threshold*/ = 100; //is_passed_theshold = TRUE
            mParameters[8] /*is_first_cycle_passed*/ = 100; //first_cycle = TRUE
        }

        if (mParameters[2] /*is_passed_threshold*/ < 50)
        {
            mParameters[2] /*is_passed_threshold*/ = 0; //is_passed_threshold = FALSE
            mParameters[11] /*is_negative_stroke_passed*/ = 0;
            if (mParameters[8] /*is_first_cycle_passed*/ > 50)
            {
                mParameters[5] /*inactive_duration*/ += ode_time_step;
            }
            /* Cheesy way of resetting the state variable */
            rDY[0] = 0;
            rDY[1] = ( 1.06669634258941e-5 - rY[1])/ode_time_step;
            rDY[2] = ( 9.57461476818094e-6 - rY[2])/ode_time_step;
            rDY[3] = ( 0.982524908182641 - rY[3])/ode_time_step;
            rDY[4] = ( 0.999999999912651 - rY[4])/ode_time_step;
            rDY[5] = ( 0.00111781215370465 - rY[5])/ode_time_step;
            rDY[6] = ( 0.769965220525522 - rY[6])/ode_time_step;
            rDY[7] = ( 0.000362322744959391 - rY[7])/ode_time_step;
            rDY[8] = ( 0.199999791554489 - rY[8])/ode_time_step;
            rDY[9] = ( 0.00441922841884097 - rY[9])/ode_time_step;
            rDY[10] = ( 0.996587027588045 - rY[10])/ode_time_step;
            rDY[11] = ( 0.0161886058733976 - rY[11])/ode_time_step;
            rDY[12] = ( 0.166158939491521 - rY[12])/ode_time_step;
            rDY[13] = ( 0.00167688179762435 - rY[13])/ode_time_step;
            rDY[14] = ( 7.61995788685097e-5 - rY[14])/ode_time_step;
            rDY[15] = ( 0.00023547935935691 - rY[15])/ode_time_step;
            rDY[16] = ( 0.00303105752591193 - rY[16])/ode_time_step;
            rDY[17] = ( 2.60090006966967 - rY[17])/ode_time_step;
            rDY[18] = ( 0.00772838990907673 - rY[18])/ode_time_step;
            rDY[19] = ( 0.101813285048293 - rY[19])/ode_time_step;
            rDY[20] = ( 0.939672327086772 - rY[20])/ode_time_step;
            rDY[21] = (163.998863583137 - rY[21])/ode_time_step;
            return;
        }

        mParameters[5] /*inactive_duration*/ = 0;

		/******************************/
		/* END */
		/******************************/

        // Inputs:
        // Time units: millisecond
        double var_chaste_interface__ICC_Membrane__Vm = (mSetVoltageDerivativeToZero ? this->mFixedVoltage : rY[0]);
        // Units: millivolt; Initial value: -67
        double var_chaste_interface__ICC_Membrane__Ca_i = rY[1];
        // Units: millimolar; Initial value: 1.06669634258941e-5
        double var_chaste_interface__d_Ltype__d_Ltype = rY[2];
        // Units: dimensionless; Initial value: 9.57461476818094e-6
        double var_chaste_interface__f_Ltype__f_Ltype = rY[3];
        // Units: dimensionless; Initial value: 0.982524908182641
        double var_chaste_interface__f_ca_Ltype__f_ca_Ltype = rY[4];
        // Units: dimensionless; Initial value: 0.999999999912651
        double var_chaste_interface__d_VDDR__d_VDDR = rY[5];
        // Units: dimensionless; Initial value: 0.00111781215370465
        double var_chaste_interface__f_VDDR__f_VDDR = rY[6];
        // Units: dimensionless; Initial value: 0.769965220525522
        double var_chaste_interface__d_CaCl__d_CaCl = rY[7];
        // Units: dimensionless; Initial value: 0.000362322744959391
        double var_chaste_interface__d_ERG__d_ERG = rY[8];
        // Units: dimensionless; Initial value: 0.199999791554489
        double var_chaste_interface__d_kv11__d_kv11 = rY[9];
        // Units: dimensionless; Initial value: 0.00441922841884097
        double var_chaste_interface__f_kv11__f_kv11 = rY[10];
        // Units: dimensionless; Initial value: 0.996587027588045
        double var_chaste_interface__d_Na__d_Na = rY[11];
        // Units: dimensionless; Initial value: 0.0161886058733976
        double var_chaste_interface__f_Na__f_Na = rY[12];
        // Units: dimensionless; Initial value: 0.166158939491521
        double var_chaste_interface__d_NSCC__d_NSCC = rY[13];
        // Units: dimensionless; Initial value: 0.00167688179762435
        double var_chaste_interface__PU_unit__Ca_PU = rY[14];
        // Units: millimolar; Initial value: 7.61995788685097e-5
        double var_chaste_interface__PU_unit__Ca_m = rY[15];
        // Units: millimolar; Initial value: 0.00023547935935691
        double var_chaste_interface__PU_unit__Ca_ER = rY[16];
        // Units: millimolar; Initial value: 0.00303105752591193
        double var_chaste_interface__PU_unit__ADP_m = rY[17];
        // Units: millimolar; Initial value: 2.60090006966967
        double var_chaste_interface__PU_unit__ADP_i = rY[18];
        // Units: millimolar; Initial value: 0.00772838990907673
        double var_chaste_interface__PU_unit__NADH_m = rY[19];
        // Units: millimolar; Initial value: 0.101813285048293
        double var_chaste_interface__PU_unit__h = rY[20];
        // Units: dimensionless; Initial value: 0.939672327086772
        double var_chaste_interface__PU_unit__deltaPsi = rY[21];
        // Units: voltage_units; Initial value: 163.998863583137


        // Mathematics
        double d_dt_chaste_interface__ICC_Membrane__Vm;
        const double var_PU_unit__J_leak = 0.01 * (var_chaste_interface__PU_unit__Ca_PU - var_chaste_interface__ICC_Membrane__Ca_i); // millimolar_per_second
        const double var_I_Ltype__I_Ltype = 2.0 * var_chaste_interface__f_Ltype__f_Ltype * var_chaste_interface__d_Ltype__d_Ltype * var_chaste_interface__f_ca_Ltype__f_ca_Ltype * (var_chaste_interface__ICC_Membrane__Vm - (log(2.5 / var_chaste_interface__ICC_Membrane__Ca_i) * 13.3568673135)); // current_units
        const double var_I_VDDR__I_VDDR = 3.0 * var_chaste_interface__f_VDDR__f_VDDR * var_chaste_interface__d_VDDR__d_VDDR * (var_chaste_interface__ICC_Membrane__Vm - (log(2.5 / var_chaste_interface__ICC_Membrane__Ca_i) * 13.3568673135)); // current_units
        const double var_J_PMCA__J_PMCA = 0.088464 / (1.0 + (0.000298 / var_chaste_interface__ICC_Membrane__Ca_i)); // millimolar_per_second
        const double var_PU_unit__J_NaCa = (0.05 * exp((var_chaste_interface__PU_unit__deltaPsi - 91.0) * 0.0187169636511)) / (1.09817777778 * (1.0 + (0.003 / var_chaste_interface__PU_unit__Ca_m))); // millimolar_per_second
        const double var_PU_unit__J_uni = ((((((0.001 * var_chaste_interface__PU_unit__Ca_PU) * 166.666666667) * pow(1.0 + (var_chaste_interface__PU_unit__Ca_PU * 166.666666667), 3.0)) / (pow(1.0 + (var_chaste_interface__PU_unit__Ca_PU * 166.666666667), 4.0) + (50.0 / pow(1.0 + (var_chaste_interface__PU_unit__Ca_PU * 2631.57894737), 2.8)))) - (var_chaste_interface__PU_unit__Ca_m * exp((var_chaste_interface__PU_unit__deltaPsi - 91.0) *  -0.0748678546044))) * (var_chaste_interface__PU_unit__deltaPsi - 91.0) * 374.339273022) / (1.0 - exp((var_chaste_interface__PU_unit__deltaPsi - 91.0) *  -0.0748678546044)); // millimolar_per_second
        const double var_PU_unit__J_SERCA = (1.8333 * pow(var_chaste_interface__PU_unit__Ca_PU, 2.0)) / (1.764e-07 + pow(var_chaste_interface__PU_unit__Ca_PU, 2.0)); // millimolar_per_second
        const double var_PU_unit__J_ERout = ((50000.0 * pow(mParameters[0] / (mParameters[0] + 0.00025), 3.0) * pow(var_chaste_interface__PU_unit__Ca_PU / (var_chaste_interface__PU_unit__Ca_PU + 0.001), 3.0) * pow(var_chaste_interface__PU_unit__h, 3.0)) + 1.666667) * (var_chaste_interface__PU_unit__Ca_ER - var_chaste_interface__PU_unit__Ca_PU); // millimolar_per_second
        const double var_PU_unit__ADP_ifree = 0.3 * var_chaste_interface__PU_unit__ADP_i; // millimolar
        const double var_PU_unit__ADP3_i = 0.45 * var_PU_unit__ADP_ifree; // millimolar
        const double var_PU_unit__ADP_mfree = 0.8 * var_chaste_interface__PU_unit__ADP_m; // millimolar
        const double var_PU_unit__ADP3_m = 0.45 * var_PU_unit__ADP_mfree; // millimolar
        const double var_PU_unit__ATP_i = 2.0 - var_chaste_interface__PU_unit__ADP_i; // millimolar
        const double var_PU_unit__ATP4_i = 0.05 * var_PU_unit__ATP_i; // millimolar
        const double var_PU_unit__ATP_m = 12.0 - var_chaste_interface__PU_unit__ADP_m; // millimolar
        const double var_PU_unit__ATP4_m = 0.05 * var_PU_unit__ATP_m; // millimolar
        const double var_PU_unit__J_ANT = (15.0 * (1.0 - (((var_PU_unit__ATP4_i * var_PU_unit__ADP3_m) / (var_PU_unit__ADP3_i * var_PU_unit__ATP4_m)) * exp(var_chaste_interface__PU_unit__deltaPsi *  -0.0374339273022)))) / ((1.0 + ((var_PU_unit__ATP4_i / var_PU_unit__ADP3_i) * exp(var_chaste_interface__PU_unit__deltaPsi *  -0.0187169636511))) * (1.0 + (var_PU_unit__ADP3_m / var_PU_unit__ATP4_m))); // millimolar_per_second
        const double var_PU_unit__f_PDHa = 1.0 / (1.0 + (1.1 * (1.0 + (15.0 / pow(1.0 + (var_chaste_interface__PU_unit__Ca_m * 20000.0), 2.0))))); // dimensionless
        const double var_PU_unit__J_glyTotal = (var_PU_unit__ATP_i * 0.13611087) / (1.0 + (4.0 * var_PU_unit__ATP_i) + ((1.0 + (2.83 * var_PU_unit__ATP_i)) * 1.3) + ((1.0 + (2.66 * var_PU_unit__ATP_i)) * 0.16)); // millimolar_per_second
        const double var_PU_unit__A_F1 = 26.7137346271 * log((1710000000.0 * var_PU_unit__ATP_m) / (var_PU_unit__ADP_mfree * 20.0)); // voltage_units
        const double var_PU_unit__A_res = 26.7137346271 * log((1.35e+18 * sqrt(var_chaste_interface__PU_unit__NADH_m)) / sqrt(8.0 - var_chaste_interface__PU_unit__NADH_m)); // voltage_units
        const double d_dt_chaste_interface__ICC_Membrane__Ca_i = 0.001 * (0.01 * (((( -1.0 * var_I_Ltype__I_Ltype) + ( -1.0 * var_I_VDDR__I_VDDR)) * 0.00740310592867) + var_PU_unit__J_leak + ( -1.0 * var_J_PMCA__J_PMCA))); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__d_Ltype__d_Ltype = 0.001 * (((1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 17.0) *  -0.232558139535))) - var_chaste_interface__d_Ltype__d_Ltype) * 381.166668969); // 'per millisecond'
        const double d_dt_chaste_interface__f_Ltype__f_Ltype = 0.001 * (((1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 43.0) * 0.112359550562))) - var_chaste_interface__f_Ltype__f_Ltype) * 4.4321705694); // 'per millisecond'
        const double d_dt_chaste_interface__f_ca_Ltype__f_ca_Ltype = 0.001 * (((1.0 - (1.0 / (1.0 + exp(((var_chaste_interface__ICC_Membrane__Ca_i - 0.0001) - 0.000214) *  -76335.8778626)))) - var_chaste_interface__f_ca_Ltype__f_ca_Ltype) * 190.583334484); // 'per millisecond'
        const double d_dt_chaste_interface__d_VDDR__d_VDDR = 0.001 * (((1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 26.0) *  -0.166666666667))) - var_chaste_interface__d_VDDR__d_VDDR) * 63.5277781614); // 'per millisecond'
        const double d_dt_chaste_interface__f_VDDR__f_VDDR = 0.001 * (((1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 66.0) * 0.166666666667))) - var_chaste_interface__f_VDDR__f_VDDR) * 9.52916672422); // 'per millisecond'
        const double d_dt_chaste_interface__d_CaCl__d_CaCl = 0.001 * (((1.0 / (1.0 + pow(0.00014 / var_chaste_interface__ICC_Membrane__Ca_i, 3.0))) - var_chaste_interface__d_CaCl__d_CaCl) * 33.3333333333); // 'per millisecond'
        const double d_dt_chaste_interface__d_ERG__d_ERG = 0.001 * (((0.2 + (0.8 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 20.0) *  -0.555555555556)))) - var_chaste_interface__d_ERG__d_ERG) * 196.770554066); // 'per millisecond'
        const double d_dt_chaste_interface__d_kv11__d_kv11 = 0.001 * (((1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 25.0) *  -0.12987012987))) - var_chaste_interface__d_kv11__d_kv11) * 118.062332439); // 'per millisecond'
        const double d_dt_chaste_interface__f_kv11__f_kv11 = 0.001 * (((0.5 + (0.5 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 44.8) * 0.227272727273)))) - var_chaste_interface__f_kv11__f_kv11) * 118.062332439); // 'per millisecond'
        const double d_dt_chaste_interface__d_Na__d_Na = 0.001 * (((1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 47.0) *  -0.208333333333))) - var_chaste_interface__d_Na__d_Na) * 103.983117932); // 'per millisecond'
        const double d_dt_chaste_interface__f_Na__f_Na = 0.001 * (((1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm + 78.0) * 0.142857142857))) - var_chaste_interface__f_Na__f_Na) * 194.968346123); // 'per millisecond'
        const double d_dt_chaste_interface__d_NSCC__d_NSCC = 0.001 * (((1.0 / (1.0 + pow(7.45e-05 / var_chaste_interface__PU_unit__Ca_PU,  -85.0))) - var_chaste_interface__d_NSCC__d_NSCC) * 2.85714285714); // 'per millisecond'
        const double d_dt_chaste_interface__PU_unit__Ca_PU = 0.001 * (0.01 * ((((var_PU_unit__J_NaCa - var_PU_unit__J_uni) * 1.2871e-13) * 1e+15) + (((var_PU_unit__J_ERout - var_PU_unit__J_SERCA) * 1e-13) * 1e+15) + ((var_PU_unit__J_leak *  -7e-13) * 1e+15))); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__PU_unit__Ca_m = 0.001 * (0.0003 * (var_PU_unit__J_uni - var_PU_unit__J_NaCa)); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__PU_unit__Ca_ER = 0.001 * (0.01 * (var_PU_unit__J_SERCA - var_PU_unit__J_ERout)); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__PU_unit__ADP_m = 0.001 * (var_PU_unit__J_ANT + ( -1.0 * (0.1111 + (0.84 * var_PU_unit__f_PDHa * var_PU_unit__J_glyTotal))) + ( -1.0 * ((((1.04489185811e-06 * exp(0.0374339273022 * var_PU_unit__A_F1)) + (exp(var_chaste_interface__PU_unit__deltaPsi * 0.112301781907) *  -1.656e-05) + (4.845e-19 * exp(0.0374339273022 * var_PU_unit__A_F1) * exp(var_chaste_interface__PU_unit__deltaPsi * 0.112301781907))) *  -0.7) / (((1.0 + (1.346e-08 * exp(0.0374339273022 * var_PU_unit__A_F1))) * 274.537838144) + ((7.739e-07 + (6.65e-15 * exp(0.0374339273022 * var_PU_unit__A_F1))) * exp(var_chaste_interface__PU_unit__deltaPsi * 0.112301781907)))))); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__PU_unit__ADP_i = 0.001 * (((var_PU_unit__J_ANT *  -1.2871e-13) * 1.42857142857e+12) + ((0.05125 * var_PU_unit__ATP_i) + 0.000109022036439) + ( -1.0 * (0.15 * var_PU_unit__J_glyTotal))); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__PU_unit__NADH_m = 0.001 * ((0.3333 + (6.3944 * var_PU_unit__f_PDHa * var_PU_unit__J_glyTotal)) - ((((2.56551579265e-12 * exp(var_PU_unit__A_res * 0.0374339273022)) + (exp(var_chaste_interface__PU_unit__deltaPsi * 0.190913029241) *  -6.394e-10) + (8.632e-27 * exp(0.0374339273022 * var_PU_unit__A_res) * exp(var_chaste_interface__PU_unit__deltaPsi * 0.190913029241))) * 0.2) / (((1.0 + (2.077e-18 * exp(0.0374339273022 * var_PU_unit__A_res))) * 75371.024573) + ((1.728e-09 + (1.059e-26 * exp(0.0374339273022 * var_PU_unit__A_res))) * exp(var_chaste_interface__PU_unit__deltaPsi * 0.190913029241))))); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__PU_unit__h = 0.001 * ((1.0 * (0.0014 - (var_chaste_interface__PU_unit__h * (var_chaste_interface__PU_unit__Ca_PU + 0.0014)))) * 0.25); // 'per millisecond'
        const double d_dt_chaste_interface__PU_unit__deltaPsi = 0.001 * ( -0.00177534422673 * ((0.0033333 * (var_chaste_interface__PU_unit__deltaPsi -  -24.6086923385)) + ( -1.0 * ((((exp(0.0374339273022 * var_PU_unit__A_res) * 2.54549724852e-12) + 7.01464834515e-16 + (exp(var_chaste_interface__PU_unit__deltaPsi * 0.190913029241) *  -6.395762e-10)) * 1.5864) / (((1.0 + (2.077e-18 * exp(0.0374339273022 * var_PU_unit__A_res))) * 75371.024573) + ((1.728e-09 + (1.059e-26 * exp(0.0374339273022 * var_PU_unit__A_res))) * exp(var_chaste_interface__PU_unit__deltaPsi * 0.190913029241))))) + var_PU_unit__J_ANT + ((((exp(0.0374339273022 * var_PU_unit__A_F1) * 1.04486536246e-06) + 2.12821912293e-08 + (exp(var_chaste_interface__PU_unit__deltaPsi * 0.112301781907) *  -1.68973e-05)) *  -2.1) / (((1.0 + (1.346e-08 * exp(0.0374339273022 * var_PU_unit__A_F1))) * 274.537838144) + ((7.739e-07 + (6.65e-15 * exp(0.0374339273022 * var_PU_unit__A_F1))) * exp(var_chaste_interface__PU_unit__deltaPsi * 0.112301781907)))) + (2.0 * var_PU_unit__J_uni))); // 'millivolt per millisecond'

        if (mSetVoltageDerivativeToZero)
        {
            d_dt_chaste_interface__ICC_Membrane__Vm = 0.0;
        }
        else
        {
            var_ICC_Membrane__Cm = 0.025; // capacitance_units
            const double var_I_Na__I_Na = 20.0 * var_chaste_interface__f_Na__f_Na * var_chaste_interface__d_Na__d_Na * (var_chaste_interface__ICC_Membrane__Vm - 40.5723805548); // current_units
            const double var_I_kv11__I_kv11 = 6.3 * var_chaste_interface__f_kv11__f_kv11 * var_chaste_interface__d_kv11__d_kv11 * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
            const double var_I_BK__I_BK = 37.3 * (1.0 / (1.0 + exp((var_chaste_interface__ICC_Membrane__Vm *  -0.0588235294118) - (2.0 * log(var_chaste_interface__ICC_Membrane__Ca_i * 1000.0))))) * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
            const double var_I_ERG__I_ERG = 2.5 * var_chaste_interface__d_ERG__d_ERG * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
            const double var_I_CaCl__I_CaCl = 10.1 * var_chaste_interface__d_CaCl__d_CaCl * (var_chaste_interface__ICC_Membrane__Vm -  -11.2332051638); // current_units
            const double var_I_NSCC__I_NSCC = 12.15 * var_chaste_interface__d_NSCC__d_NSCC * (var_chaste_interface__ICC_Membrane__Vm - 4.40291009976e-06); // current_units
            const double var_I_bk__I_bk = 0.15 * (var_chaste_interface__ICC_Membrane__Vm -  -75.909256616); // current_units
            d_dt_chaste_interface__ICC_Membrane__Vm = 0.001 * (( -1.0 / var_ICC_Membrane__Cm) * (var_I_Na__I_Na + var_I_Ltype__I_Ltype + var_I_VDDR__I_VDDR + var_I_kv11__I_kv11 + var_I_ERG__I_ERG + var_I_BK__I_BK + var_I_CaCl__I_CaCl + var_I_NSCC__I_NSCC + var_I_bk__I_bk + (var_J_PMCA__J_PMCA * 135.07844))); // 'millivolt per millisecond'
        }

        rDY[0] = d_dt_chaste_interface__ICC_Membrane__Vm;
        rDY[1] = d_dt_chaste_interface__ICC_Membrane__Ca_i;
        rDY[2] = d_dt_chaste_interface__d_Ltype__d_Ltype;
        rDY[3] = d_dt_chaste_interface__f_Ltype__f_Ltype;
        rDY[4] = d_dt_chaste_interface__f_ca_Ltype__f_ca_Ltype;
        rDY[5] = d_dt_chaste_interface__d_VDDR__d_VDDR;
        rDY[6] = d_dt_chaste_interface__f_VDDR__f_VDDR;
        rDY[7] = d_dt_chaste_interface__d_CaCl__d_CaCl;
        rDY[8] = d_dt_chaste_interface__d_ERG__d_ERG;
        rDY[9] = d_dt_chaste_interface__d_kv11__d_kv11;
        rDY[10] = d_dt_chaste_interface__f_kv11__f_kv11;
        rDY[11] = d_dt_chaste_interface__d_Na__d_Na;
        rDY[12] = d_dt_chaste_interface__f_Na__f_Na;
        rDY[13] = d_dt_chaste_interface__d_NSCC__d_NSCC;
        rDY[14] = d_dt_chaste_interface__PU_unit__Ca_PU;
        rDY[15] = d_dt_chaste_interface__PU_unit__Ca_m;
        rDY[16] = d_dt_chaste_interface__PU_unit__Ca_ER;
        rDY[17] = d_dt_chaste_interface__PU_unit__ADP_m;
        rDY[18] = d_dt_chaste_interface__PU_unit__ADP_i;
        rDY[19] = d_dt_chaste_interface__PU_unit__NADH_m;
        rDY[20] = d_dt_chaste_interface__PU_unit__h;
        rDY[21] = d_dt_chaste_interface__PU_unit__deltaPsi;

		/******************************/
		/* FSM modification Ca2+ Tracker */
		/******************************/

        if (rY[0] > -30)
	    {
            mParameters[11] /*is_negative_stroke_passed*/ = 100; // Peak tracked
            return;
        }

        if (mParameters[11] /*is_negative_stroke_passed*/ > 50 && ((rY[1] + rDY[1]*mParameters[7] /*ode-time-step*/) <  1.26e-5)) /* Negative stroke passed and [Ca2+]i is 0.*/
        {
            mParameters[2] /*is_passed_threshold*/ = 0; //is_passed_threshold = FALSE
            mParameters[11] /*is_negative_stroke_passed*/ = 0;
            if (mParameters[8] /*is_first_cycle_passed*/ > 50)
		    {
                mParameters[5] /*inactive_duration*/ += ode_time_step;
		    }

            /* Cheesy way of resetting the state variable */
            rDY[0] = 0;
            rDY[1] = ( 1.06669634258941e-5 - rY[1])/ode_time_step;
            rDY[2] = ( 9.57461476818094e-6 - rY[2])/ode_time_step;
            rDY[3] = ( 0.982524908182641 - rY[3])/ode_time_step;
            rDY[4] = ( 0.999999999912651 - rY[4])/ode_time_step;
            rDY[5] = ( 0.00111781215370465 - rY[5])/ode_time_step;
            rDY[6] = ( 0.769965220525522 - rY[6])/ode_time_step;
            rDY[7] = ( 0.000362322744959391 - rY[7])/ode_time_step;
            rDY[8] = ( 0.199999791554489 - rY[8])/ode_time_step;
            rDY[9] = ( 0.00441922841884097 - rY[9])/ode_time_step;
            rDY[10] = ( 0.996587027588045 - rY[10])/ode_time_step;
            rDY[11] = ( 0.0161886058733976 - rY[11])/ode_time_step;
            rDY[12] = ( 0.166158939491521 - rY[12])/ode_time_step;
            rDY[13] = ( 0.00167688179762435 - rY[13])/ode_time_step;
            rDY[14] = ( 7.61995788685097e-5 - rY[14])/ode_time_step;
            rDY[15] = ( 0.00023547935935691 - rY[15])/ode_time_step;
            rDY[16] = ( 0.00303105752591193 - rY[16])/ode_time_step;
            rDY[17] = ( 2.60090006966967 - rY[17])/ode_time_step;
            rDY[18] = ( 0.00772838990907673 - rY[18])/ode_time_step;
            rDY[19] = ( 0.101813285048293 - rY[19])/ode_time_step;
            rDY[20] = ( 0.939672327086772 - rY[20])/ode_time_step;
            rDY[21] = (163.998863583137 - rY[21])/ode_time_step;
            return;
        }
		/******************************/
		/* END */
		/******************************/
    }


	template<>
	void OdeSystemInformation<CellICCBioPhy>::Initialise(void)
	{
    	this->mSystemName = "ICC_model";
    	this->mFreeVariableName = "Time__time";
    	this->mFreeVariableUnits = "millisecond";

    	this->mVariableNames.push_back("ICC_Membrane__Vm");
    	this->mVariableUnits.push_back("millivolt");
    	this->mInitialConditions.push_back(-67);

    	this->mVariableNames.push_back("ICC_Membrane__Ca_i");
    	this->mVariableUnits.push_back("millimolar");
    	this->mInitialConditions.push_back(1.06669634258941e-5);

    	this->mVariableNames.push_back("d_Ltype__d_Ltype");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(9.57461476818094e-6);

    	this->mVariableNames.push_back("f_Ltype__f_Ltype");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.982524908182641);

    	this->mVariableNames.push_back("f_ca_Ltype__f_ca_Ltype");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.999999999912651);

    	this->mVariableNames.push_back("d_VDDR__d_VDDR");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.00111781215370465);

    	this->mVariableNames.push_back("f_VDDR__f_VDDR");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.769965220525522);

    	this->mVariableNames.push_back("d_CaCl__d_CaCl");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.000362322744959391);

    	this->mVariableNames.push_back("d_ERG__d_ERG");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.199999791554489);

    	this->mVariableNames.push_back("d_kv11__d_kv11");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.00441922841884097);

    	this->mVariableNames.push_back("f_kv11__f_kv11");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.996587027588045);

    	this->mVariableNames.push_back("d_Na__d_Na");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.0161886058733976);

    	this->mVariableNames.push_back("f_Na__f_Na");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.166158939491521);

    	this->mVariableNames.push_back("d_NSCC__d_NSCC");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.00167688179762435);

    	this->mVariableNames.push_back("PU_unit__Ca_PU");
    	this->mVariableUnits.push_back("millimolar");
    	this->mInitialConditions.push_back(7.61995788685097e-5);

    	this->mVariableNames.push_back("PU_unit__Ca_m");
    	this->mVariableUnits.push_back("millimolar");
    	this->mInitialConditions.push_back(0.00023547935935691);

    	this->mVariableNames.push_back("PU_unit__Ca_ER");
    	this->mVariableUnits.push_back("millimolar");
    	this->mInitialConditions.push_back(0.00303105752591193);

    	this->mVariableNames.push_back("PU_unit__ADP_m");
    	this->mVariableUnits.push_back("millimolar");
    	this->mInitialConditions.push_back(2.60090006966967);

    	this->mVariableNames.push_back("PU_unit__ADP_i");
    	this->mVariableUnits.push_back("millimolar");
    	this->mInitialConditions.push_back(0.00772838990907673);

    	this->mVariableNames.push_back("PU_unit__NADH_m");
	  	this->mVariableUnits.push_back("millimolar");
    	this->mInitialConditions.push_back(0.101813285048293);

    	this->mVariableNames.push_back("PU_unit__h");
    	this->mVariableUnits.push_back("dimensionless");
    	this->mInitialConditions.push_back(0.939672327086772);

    	this->mVariableNames.push_back("PU_unit__deltaPsi");
    	this->mVariableUnits.push_back("voltage_units");
    	this->mInitialConditions.push_back(163.998863583137);

    	this->mParameterNames.push_back("IP3Par");
    	this->mParameterUnits.push_back("millimolar");

    	/******************************/
		/* FSM modification parameters */
		/******************************/
    	this->mParameterNames.push_back("V_excitation");
    	this->mParameterUnits.push_back("millivolt");

    	this->mParameterNames.push_back("is_passed_threshold");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("active_duration");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("dead_time");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("inactive_duration");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("live_time");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("ode_time_step");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("first_cycle_passed");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("t_start");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("t_elapsed");
    	this->mParameterUnits.push_back("dimensionless");

    	this->mParameterNames.push_back("negative_stroke_passed");
    	this->mParameterUnits.push_back("dimensionless");
    	/******************************/
		/* END */
		/******************************/

    	this->mInitialised = true;
	}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CellICCBioPhy)

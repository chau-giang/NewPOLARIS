#include "Grid.h"
#include <CCfits/CCfits>
#include <cmath>

void CGridBasic::updateDataRange(cell_basic * cell, parameters & param)
{
    double gas_dens = 0;
    double dust_dens = 0;
    double mx = 0;
    double my = 0;
    double mz = 0;
    double vx = 0;
    double vy = 0;
    double vz = 0;
    double px = 0;
    double py = 0;
    double pz = 0;
    double a_alg = 0;
    double amaxJB_Lar = 0;
    double akrat_lowJ = 0;
    double akrat_highJ = 0;
    double adisr = 0;
    double max_adisr = 0;
    double size_param_modif = 0;
    double abar_low_lower = 0;
    double abar_low_upper = 0; 
    double abar_high_lower = 0;
    double abar_high_upper = 0;
    double adg_lower = 0;
    double adg_upper = 0;
    double adg_10_lower = 0;
    double adg_10_upper = 0;
    double abs_ini_temp = 0;
    double dust_temp = 0;
    double abs_ini = 0;
    double gas_temp = 0;
    double mach = 0;
    double delta = 0;
     

    if(!data_pos_gd_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_gd_list.size(); i_dens++)
            cell->convertData(data_pos_gd_list[i_dens], conv_dens_in_SI);
        gas_dens = getGasDensity(cell);

        if(gas_dens > max_gas_dens)
            max_gas_dens = gas_dens;
        if(gas_dens < min_gas_dens)
            min_gas_dens = gas_dens;
    }

    if(!data_pos_dd_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_dd_list.size(); i_dens++)
            cell->convertData(data_pos_dd_list[i_dens], conv_dens_in_SI);

        dust_dens = getDustDensity(cell);

        if(dust_dens > max_dust_dens)
            max_dust_dens = dust_dens;
        if(dust_dens < min_dust_dens)
            min_dust_dens = dust_dens;
    }

    if(!data_pos_dt_list.empty())
    {
        dust_temp = cell->getData(data_pos_dt_list[0]);
        // to do if conversion is implemented

        if(dust_temp > max_dust_temp)
            max_dust_temp = dust_temp;
        if(dust_temp < min_dust_temp)
            min_dust_temp = dust_temp;
    }

    if(data_pos_tg != MAX_UINT)
    {
        gas_temp = cell->getData(data_pos_tg);
        // to do if conversion is implemented
        if(gas_temp > max_gas_temp)
            max_gas_temp = gas_temp;
        if(gas_temp < min_gas_temp)
            min_gas_temp = gas_temp;
    }

    if(data_pos_mx != MAX_UINT)
    {
        cell->convertData(data_pos_mx, conv_Bfield_in_SI);
        mx = cell->getData(data_pos_mx);
    }

    if(data_pos_my != MAX_UINT)
    {
        cell->convertData(data_pos_my, conv_Bfield_in_SI);
        my = cell->getData(data_pos_my);
    }

    if(data_pos_mz != MAX_UINT)
    {
        cell->convertData(data_pos_mz, conv_Bfield_in_SI);
        mz = cell->getData(data_pos_mz);
    }

    if(data_pos_vx != MAX_UINT)
    {
        cell->convertData(data_pos_vx, conv_Vfield_in_SI);
        vx = cell->getData(data_pos_vx);
    }

    if(data_pos_vy != MAX_UINT)
    {
        cell->convertData(data_pos_vy, conv_Vfield_in_SI);
        vy = cell->getData(data_pos_vy);
    }

    if(data_pos_vz != MAX_UINT)
    {
        cell->convertData(data_pos_vz, conv_Vfield_in_SI);
        vz = cell->getData(data_pos_vz);
    }

    if(!data_pos_aalg_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_aalg_list.size(); i_dens++)
        {
            a_alg = cell->getData(data_pos_aalg_list[i_dens]);

            if(a_alg > float(aalg_max))
                aalg_max = (double)a_alg;

            if(a_alg < float(aalg_min))
                aalg_min = (double)a_alg;
        }
    }

    if(data_pos_amin != MAX_UINT)
    {
        double a_min = cell->getData(data_pos_amin);

        if(a_min > float(a_min_max))
            a_min_max = a_min;

        if(a_min < float(a_min_min))
            a_min_min = a_min;
    }

    if(data_pos_amax != MAX_UINT)
    {
        double a_max = cell->getData(data_pos_amax);

        if(a_max > float(a_max_max))
            a_max_max = a_max;

        if(a_max < float(a_max_min))
            a_max_min = a_max;
    }

    if(data_pos_size_param != MAX_UINT)
    {
        uint size_param = cell->getData(data_pos_size_param);

        if(size_param > float(size_param_max))
            size_param_max = size_param;

        if(size_param < float(size_param_min))
            size_param_min = size_param;
    }

    if(data_pos_id != MAX_UINT)
    {
        uint dust_id = cell->getData(data_pos_id);

        if(dust_id > float(dust_id_max))
            dust_id_max = dust_id;

        if(dust_id < float(dust_id_min))
            dust_id_min = dust_id;
    }

    // data positions for synchrotron
    if(data_pos_n_th != MAX_UINT)
    {
        cell->convertData(data_pos_n_th, conv_dens_in_SI);
        double data = cell->getData(data_pos_n_th);

        if(data < min_n_th)
            min_n_th = data;

        if(data > max_n_th)
            max_n_th = data;
    }

    if(data_pos_T_e != MAX_UINT)
    {
        double data = cell->getData(data_pos_T_e);

        if(data < min_T_e)
            min_T_e = data;

        if(data > max_T_e)
            max_T_e = data;
    }

    if(data_pos_n_cr != MAX_UINT)
    {
        cell->convertData(data_pos_n_cr, conv_dens_in_SI);
        double data = cell->getData(data_pos_n_cr);

        if(data < min_n_cr)
            min_n_cr = data;

        if(data > max_n_cr)
            max_n_cr = data;
    }

    if(data_pos_g_min != MAX_UINT)
    {
        double data = cell->getData(data_pos_g_min);

        if(data < min_g_min)
            min_g_min = data;

        if(data > max_g_min)
            max_g_min = data;
    }

    if(data_pos_g_max != MAX_UINT)
    {
        double data = cell->getData(data_pos_g_max);

        if(data < min_g_max)
            min_g_max = data;

        if(data > max_g_max)
            max_g_max = data;
    }

    if(data_pos_p != MAX_UINT)
    {
        double data = cell->getData(data_pos_p);

        if(data < min_p)
            min_p = data;

        if(data > max_p)
            max_p = data;
    }


    if(!data_pos_adisr_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_adisr_list.size(); i_dens++)
        {
            adisr = cell->getData(data_pos_adisr_list[i_dens]);

            if(adisr > float(adisr_max))
                adisr_max = (double)adisr;

            if(adisr < float(adisr_min))
                adisr_min = (double)adisr;
        }
    }

    if(!data_pos_max_adisr_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_max_adisr_list.size(); i_dens++)
        {
            max_adisr = cell->getData(data_pos_max_adisr_list[i_dens]);

            if(max_adisr > float(max_adisr_max))
                max_adisr_max = (double)max_adisr;

            if(max_adisr < float(max_adisr_min))
                max_adisr_min = (double)max_adisr;
        }
    }

    if(!data_pos_param_modif_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_param_modif_list.size(); i_dens++)
        {
            size_param_modif = cell->getData(data_pos_param_modif_list[i_dens]);

            if(size_param_modif > float(size_param_modif_max))
                size_param_modif_max = (double)size_param_modif;

            if(size_param_modif < float(size_param_modif_min))
                size_param_modif_min = (double)size_param_modif;
        }
    }

    if(!data_pos_barnet_low_J_lower_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_barnet_low_J_lower_list.size(); i_dens++)
        {
            abar_low_lower = cell->getData(data_pos_barnet_low_J_lower_list[i_dens]);

            if(abar_low_lower > float(abar_low_lower_max))
                abar_low_lower_max = (double)abar_low_lower;

            if(abar_low_lower < float(abar_low_lower_min))
                abar_low_lower_min = (double)abar_low_lower;
        }
    }
    
    if(!data_pos_barnet_low_J_upper_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_barnet_low_J_upper_list.size(); i_dens++)
        {
            abar_low_upper = cell->getData(data_pos_barnet_low_J_upper_list[i_dens]);

            if(abar_low_upper > float(abar_low_upper_max))
                abar_low_upper_max = (double)abar_low_upper;

            if(abar_low_upper < float(abar_low_upper_min))
                abar_low_upper_min = (double)abar_low_upper;
        }
    } 
 
    if(!data_pos_barnet_high_J_lower_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_barnet_high_J_lower_list.size(); i_dens++)
        {
            abar_high_lower = cell->getData(data_pos_barnet_high_J_lower_list[i_dens]);

            if(abar_high_lower > float(abar_high_lower_max))
                abar_high_lower_max = (double)abar_high_lower;

            if(abar_high_lower < float(abar_high_lower_min))
                abar_high_lower_min = (double)abar_high_lower;
        }
    }
    
    if(!data_pos_barnet_high_J_upper_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_barnet_high_J_upper_list.size(); i_dens++)
        {
            abar_high_upper = cell->getData(data_pos_barnet_high_J_upper_list[i_dens]);

            if(abar_high_upper > float(abar_high_upper_max))
                abar_high_upper_max = (double)abar_high_upper;

            if(abar_high_upper < float(abar_high_upper_min))
                abar_high_upper_min = (double)abar_high_upper;
        }
    }
    
    if(!data_pos_dg_lower_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_dg_lower_list.size(); i_dens++)
        {
            adg_lower = cell->getData(data_pos_dg_lower_list[i_dens]);

            if(adg_lower > float(adg_lower_max))
                adg_lower_max = (double)adg_lower;

            if(adg_lower < float(adg_lower_min))
                adg_lower_min = (double)adg_lower;
        }
    }
    
    if(!data_pos_dg_upper_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_dg_upper_list.size(); i_dens++)
        {
            adg_upper = cell->getData(data_pos_dg_upper_list[i_dens]);

            if(adg_upper > float(adg_upper_max))
                adg_upper_max = (double)adg_upper;

            if(adg_upper < float(adg_upper_min))
                adg_upper_min = (double)adg_upper;
        }
    }
    
    if(!data_pos_dg_10_lower_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_dg_10_lower_list.size(); i_dens++)
        {
            adg_10_lower = cell->getData(data_pos_dg_10_lower_list[i_dens]);

            if(adg_10_lower > float(adg_10_lower_max))
                adg_10_lower_max = (double)adg_10_lower;

            if(adg_10_lower < float(adg_10_lower_min))
                adg_10_lower_min = (double)adg_10_lower;
        }
    }
    
    if(!data_pos_dg_10_upper_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_dg_10_upper_list.size(); i_dens++)
        {
            adg_10_upper = cell->getData(data_pos_dg_10_upper_list[i_dens]);

            if(adg_10_upper > float(adg_10_upper_max))
                adg_10_upper_max = (double)adg_10_upper;

            if(adg_10_upper < float(adg_10_upper_min))
                adg_10_upper_min = (double)adg_10_upper;
        }
    }
    
    // conversion from initial dust temperature to default absorption rate in the grid
    if(!data_pos_abs_ini_list.empty())
    {
        abs_ini = cell->getData(data_pos_abs_ini_list[0]);
        // to do if conversion is implemented

        if(abs_ini > max_abs_ini)
            max_abs_ini = abs_ini;
        if(abs_ini < min_abs_ini)
            min_abs_ini = abs_ini;
    }


    if(!data_pos_amaxJB_Lar_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_amaxJB_Lar_list.size(); i_dens++)
        {
            amaxJB_Lar = cell->getData(data_pos_amaxJB_Lar_list[i_dens]);

            if(amaxJB_Lar > float(max_amaxJB_Lar))
                max_amaxJB_Lar = (double)amaxJB_Lar;

            if(amaxJB_Lar < float(min_amaxJB_Lar))
                min_amaxJB_Lar = (double)amaxJB_Lar;
        }
    }

    if(!data_pos_akrat_lowJ_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_akrat_lowJ_list.size(); i_dens++)
        {
            akrat_lowJ = cell->getData(data_pos_akrat_lowJ_list[i_dens]);

            if(akrat_lowJ > float(max_akrat_lowJ))
                max_akrat_lowJ = (double)akrat_lowJ;

            if(akrat_lowJ < float(min_akrat_lowJ))
                min_akrat_lowJ = (double)akrat_lowJ;
        }
    }

    if(!data_pos_akrat_highJ_list.empty())
    {
        for(uint i_dens = 0; i_dens < data_pos_akrat_highJ_list.size(); i_dens++)
        {
            akrat_highJ = cell->getData(data_pos_akrat_highJ_list[i_dens]);

            if(akrat_highJ > float(max_akrat_highJ))
                max_akrat_highJ = (double)akrat_highJ;

            if(akrat_highJ < float(min_akrat_highJ))
                min_akrat_highJ = (double)akrat_highJ;
        }
    }
    
 
    double Bfield = sqrt(mx * mx + my * my + mz * mz);
    double Vfield = sqrt(vx * vx + vy * vy + vz * vz);

    if(Bfield > 0)
    {
        if(dust_temp * gas_temp * gas_dens >= 0)
        {
            delta = CMathFunctions::calc_delta(Bfield, dust_temp, gas_temp, gas_dens) * delta0;
        }
    }
    else
    {
        Bfield = 0;
        delta = 0;
    }

    if(delta > max_delta)
        max_delta = delta;
    if(delta < min_delta)
        min_delta = delta;

    if(Bfield > max_mag)
        max_mag = Bfield;
    if(Bfield < min_mag)
        min_mag = Bfield;

    meanBdir += Vector3D(mx, my, mz);

    if(Vfield >= 0)
    {
        if(gas_temp > 0)
            mach = Vfield / sqrt(con_kB * gas_temp / (mu * m_H));
        else
            mach = 0;
    }

    if(Vfield > max_vel)
        max_vel = Vfield;
    if(Vfield < min_vel)
        min_vel = Vfield;
    if(mach > max_mach)
        max_mach = mach;
    if(mach < min_mach)
        min_mach = mach;

    meanVdir += Vector3D(vx, vy, vz);
}

bool CGridBasic::fillGridWithOpiateData(uint col_id)
{
    /* uint cell_count = 0;
     uint found_count = 0;
     //#pragma omp parallel for schedule(dynamic)
     for(long i_cell = 0; i_cell < long(max_cells); i_cell++)
     {
 #pragma omp critical
         {
             if(cell_count % 500 == 0)
             {
                 cout << "-> Filling grid with OPIATE data  : "
                         << 100 * float(cell_count) / float(max_cells) << " [%] \r";
             }
         }

         cell_count++;

         cell_basic * cell = cell_list[i_cell];

         uint id = getOpiateID(cell);
         double val = 0;

         if(id != MAX_UINT)
         {
             val = opiate->getData(id, col_id);
             setOpiateTestData(cell, val);
             found_count++;
         }
         else
             setOpiateTestData(cell, 0);
     }

     cout << CLR_LINE;
     cout << " - " << found_count << " of " << max_cells << " cells match with the OPIATE
 paramter file." << endl;
     */
    return true;
}

uint CGridBasic::validateDataPositions(parameters & param)
{
    uint tmp_data_offset = 0;

    cout << CLR_LINE;

    if(data_pos_gd_list.size() == 0)
    {
        cout << "\nERROR: Grid contains no gas (number) density!" << endl;
        cout << "       No RT calculations possible!" << endl;
        return MAX_UINT;
    }

    if(nr_mixtures > 0 && (param.isMonteCarloSimulation() || param.getCommand() == CMD_DUST_EMISSION ||
                           param.getCommand() == CMD_LINE_EMISSION || param.getCommand() == CMD_FORCE ||
                           param.getCommand() == CMD_PROBING))
    {
        // Get Number of density fields for temperature calculation
        if(!data_pos_dd_list.empty())
            nr_densities = data_pos_dd_list.size();
        else
            nr_densities = data_pos_gd_list.size();

        // Precalculate the number of temperature entries, if the grid has a
        // temperature for each grain size or stochastically heated grains
        for(uint i_density = 0; i_density < nr_densities; i_density++)
        {
            multi_temperature_entries += nr_dust_temp_sizes[i_density] + 1;
            stochastic_temperature_entries += nr_stochastic_sizes[i_density] + 1;
            wavelength_entries += nr_wavelength[i_density] + 1;
        }

        // Check for a valid combination between densities and dust mixtures
        if(nr_densities > 1 && nr_mixtures < nr_densities)
        {
            cout << "\nERROR: Amount of densities in the grid (" << nr_densities
                 << ") does not fit with the defined dust mixtures (" << nr_mixtures << ")!\n"
                 << "(Use a grid with only one density distribution or define more/less "
                    "dust mixtures!)"
                 << endl;
            return MAX_UINT;
        }

        // Init list to know how many dust sizes are used per dust component
        size_skip = new uint[nr_densities];

        // Calculate the entries for the temperature that have to be added
        if(param.getDustTempMulti())
            for(uint i_density = 0; i_density < nr_densities; i_density++)
                size_skip[i_density] = nr_dust_temp_sizes[i_density]; //N_temp
                
        else if(param.getStochasticHeatingMaxSize() > 0 && !param.getSaveRadiationField())
            for(uint i_density = 0; i_density < nr_densities; i_density++)
                size_skip[i_density] = nr_stochastic_sizes[i_density]; // N_a
        else
            for(uint i_density = 0; i_density < nr_densities; i_density++)
                size_skip[i_density] = 1;
 
    }

    switch(param.getCommand())
    {
        case CMD_SYNCHROTRON:
            if(CheckSynchrotron(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_OPIATE:
            if(CheckOpiate(param) == MAX_UINT)
                return MAX_UINT;
            break;

        // ---------------------------------------------------
        case CMD_TEMP:
            if(CheckTemp(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_TEMP_RAT:
            if(CheckTemp(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            if(CheckRat(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_TEMP_DISR:
            if(CheckTemp(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            if(CheckRATD(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            break;

        case CMD_TEMP_RAT_DISR:
            if(CheckTemp(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            if(CheckRat(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            if(CheckRATD(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        // ---------------------------------------------------
        case CMD_RAT:
            if(CheckRat(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_DISR:
            if(CheckRATD(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_RAT_DISR:
            if(CheckRat(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;

            if(CheckRATD(param, tmp_data_offset) == MAX_UINT)
                return MAX_UINT;
            break;

        // ---------------------------------------------------
        case CMD_DUST_EMISSION:
            if(CheckDustEmission(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_DUST_SCATTERING:
            if(CheckDustScattering(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_FORCE:
            if(CheckRadiationForce(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_LINE_EMISSION:
            if(CheckLineEmission(param) == MAX_UINT)
                return MAX_UINT;
            break;

        case CMD_PROBING:
            if(CheckProbing(param) == MAX_UINT)
                return MAX_UINT;
            break;

        default:
            cout << "\nERROR: Command is unknown!" << endl;
            return MAX_UINT;
    }

    return tmp_data_offset;
}

void CGridBasic::printPhysicalParameters()
{
    cout << "- Volume (total, cells, diff)   : " << total_volume << " [m^3], " << cell_volume << " [m^3], "
         << float(100.0 * abs(total_volume - cell_volume) / max(total_volume, cell_volume)) << " [%]" << endl;
    cout << "- Total gas mass                : " << total_gas_mass / M_sun << " [M_sun], " << total_gas_mass
         << " [kg]" << endl;
    cout << "- Grid length         (min,max) : [" << min_len << ", " << max_len << "] [m]" << endl;
    if(gas_is_mass_density)
        cout << "- Gas mass density    (min,max) : [" << min_gas_dens << ", " << max_gas_dens << "] [kg m^-3]"
             << endl;
    else
        cout << "- Gas number density  (min,max) : [" << min_gas_dens << ", " << max_gas_dens << "] [m^-3]"
             << endl;
    if(!data_pos_dd_list.empty())
        if(dust_is_mass_density)
            cout << "- Dust mass density   (min,max) : [" << min_dust_dens << ", " << max_dust_dens
                 << "] [kg m^-3]" << endl;
        else
            cout << "- Dust number density (min,max) : [" << min_dust_dens << ", " << max_dust_dens
                 << "] [m^-3]" << endl;
    if(data_pos_tg != MAX_UINT)
        cout << "- Gas temperature     (min,max) : [" << min_gas_temp << ", " << max_gas_temp << "] [K]"
             << endl;
    else
        cout << "- Gas temperature     (min,max) : none" << endl;

    if(!data_pos_dt_list.empty())
        cout << "- Dust temperature    (min,max) : [" << min_dust_temp << ", " << max_dust_temp << "] [K]"
             << endl;
    else
        cout << "- Dust temperature    (min,max) : none" << endl;

    if(data_pos_mx != MAX_UINT)
    {
        meanBdir.normalize();
        cout << "- Magnetic field      (min,max) : [" << min_mag << ", " << max_mag << "] [T]" << endl;
        cout << "- Mean direction      (norm.)   : X: " << meanBdir.X() << " Y: " << meanBdir.Y()
             << " Z: " << meanBdir.Z() << endl;
        //cout << "- Delta0              (min,max) : [" << min_delta << ", " << max_delta << "] [m]" << endl;
             
    }
    else
        cout << "- Magnetic field      (min,max) : none" << endl;

    if(!data_pos_aalg_list.empty())
        cout << "- a_alig              (min,max) : [" << aalg_min << ", " << aalg_max << "] [m]" << endl;

    if(!data_pos_amaxJB_Lar_list.empty())
        cout << "- amaxJB_Lar         (min,max) : [" << min_amaxJB_Lar << ", " << max_amaxJB_Lar << "] [m]" << endl;
        
    if(!data_pos_adisr_list.empty())
        cout << "- a_disr              (min,max) : [" << adisr_min << ", " << adisr_max << "] [m]" << endl;

    if(!data_pos_max_adisr_list.empty())
        cout << "- a_disr_max           (min,max) : [" << max_adisr_min << ", " << max_adisr_max << "] [m]" << endl;

    if(!data_pos_param_modif_list.empty())
        cout << "- size distribution modify  (min,max) : [" << size_param_modif_min << ", " << size_param_modif_max << "]" << endl;
        
    if(!data_pos_barnet_low_J_lower_list.empty())
        cout << "- amin_aJ_lowJ       (min,max) : [" << abar_low_lower_min << ", " << abar_low_lower_max << "] [m]" << endl;
        
    if(!data_pos_barnet_low_J_upper_list.empty())
        cout << "- amax_aJ_lowJ    (min,max) : [" << abar_low_upper_min << ", " << abar_low_upper_max << "] [m]" << endl;

 
    if(!data_pos_barnet_high_J_lower_list.empty())
        cout << "- amin_aJ_highJ   (min,max) : [" << abar_high_lower_min << ", " << abar_high_lower_max << "] [m]" << endl;
        
    if(!data_pos_barnet_high_J_upper_list.empty())
        cout << "- amax_aJ_highJ   (min,max) : [" << abar_high_upper_min << ", " << abar_high_upper_max << "] [m]" << endl;
        
        
    //if(!data_pos_dg_lower_list.empty())
    //    cout << "- aminJB_DG_0.5   (min,max) : [" << adg_lower_min << ", " << adg_lower_max << "] [m]" << endl;
        
    if(!data_pos_dg_upper_list.empty())
        cout << "- amaxJB_DG_0.5   (min,max) : [" << adg_upper_min << ", " << adg_upper_max << "] [m]" << endl;
        
    //if(!data_pos_dg_10_lower_list.empty())
    //    cout << "- aminJB_DG_1   (min,max) : [" << adg_10_lower_min << ", " << adg_10_lower_max << "] [m]" << endl;
        
    if(!data_pos_dg_10_upper_list.empty())
        cout << "- amaxJB_DG_1   (min,max) : [" << adg_10_upper_min << ", " << adg_10_upper_max << "] [m]" << endl;
        
    if(!data_pos_akrat_lowJ_list.empty())
        cout << "- akrat_lowJ   (min,max) : [" << min_akrat_lowJ << ", " << max_akrat_lowJ << "] [m]" << endl;
        
    if(!data_pos_akrat_highJ_list.empty())
        cout << "- akrat_highJ   (min,max) : [" << min_akrat_highJ << ", " << max_akrat_highJ << "] [m]" << endl;
        

    if(data_pos_vx != MAX_UINT)
    {
        meanVdir.normalize();
        cout << "- Velocity field      (min,max) : [" << min_vel << ", " << max_vel << "] [m/s]" << endl;
        cout << "- Mean direction      (norm.)   : X: " << meanVdir.X() << " Y: " << meanVdir.Y()
             << " Z: " << meanVdir.Z() << endl;
        cout << "- Mach number         (min,max) : [" << min_mach << ", " << max_mach << "]" << endl;
    }

    if(data_pos_amin != MAX_UINT)
        cout << "- Minimum grain size  (min,max) : [" << a_min_min << ", " << a_min_max << "] [m]" << endl;

    if(data_pos_amax != MAX_UINT)
        cout << "- Maximum grain size  (min,max) : [" << a_max_min << ", " << a_max_max << "] [m]" << endl;

    if(data_pos_size_param != MAX_UINT)
        cout << "- Dust size parameter (min,max) : [" << size_param_min << ", " << size_param_max << "]"
             << endl;

    if(data_pos_id != MAX_UINT)
        cout << "- Dust mixture ID     (min,max) : [" << dust_id_min << ", " << dust_id_max << "]" << endl;

    if(data_pos_n_cr != MAX_UINT)
    {
        cout << "- CR el. density      (min,max) : [" << min_n_cr << "; " << max_n_cr << "] [m^-3]" << endl;

        if(data_pos_g_min != MAX_UINT)
            cout << "- Gamma_min           (min,max) : [" << min_g_min << "; " << max_g_min << "]" << endl;

        if(data_pos_g_max != MAX_UINT)
            cout << "- Gamma_max           (min,max) : [" << min_g_max << "; " << max_g_max << "]" << endl;

        if(data_pos_p != MAX_UINT)
            cout << "- El. energy index p  (min,max) : [" << min_p << "; " << max_p << "]" << endl;
    }
    else
        cout << "- CR el. density      (min,max) : none   " << endl;

    if(data_pos_n_th != MAX_UINT)
    {
        cout << "- Therm. el. density  (min,max) : [" << min_n_th << "; " << max_n_th << "] [m]" << endl;

        if(data_pos_T_e != MAX_UINT)
        {
            if(min_T_e == 1e300)
                cout << "- Electron temperature          : same as dust temperature" << endl;
            else
                cout << "- Electron temp.      (min,max) : [" << min_T_e << "; " << max_T_e << "] [K]"
                     << endl;
        }
    }
    else
        cout << "- Therm. el. density  (min,max) : none" << endl;

    if(nrOfOpiateIDs > 0 || nrOfDensRatios > 0)
    {
        cout << SEP_LINE;
        cout << "Additional grid data:" << endl;
    }

    if(nrOfDensRatios > 0)
    {
        cout << "- Density. ratio IDs: ";
        cout << 1 << ":" << pos_GasSpecRatios[0];

        for(uint i = 1; i < nrOfDensRatios; i++)
            cout << ", " << i + 1 << ":" << pos_GasSpecRatios[i];

        cout << endl;
    }

    if(nrOfOpiateIDs > 0)
    {
        cout << "- Unique OPIATE IDs : ";
        cout << 1 << ":" << pos_OpiateIDS[0];

        for(uint i = 1; i < nrOfOpiateIDs; i++)
            cout << ", " << i + 1 << ":" << pos_OpiateIDS[i];

        cout << endl;
    }

    if(data_pos_op != UINT_MAX)
        cout << " - Unique OPIATE IDs" << endl;
};

bool CGridBasic::writeAMIRAFiles(string path, parameters & param, uint bins)
{
    if(bins == 0)
        return true;

    bool plt_gas_dens = (!data_pos_gd_list.empty()) && param.isInPlotList(GRIDgas_dens);
    // bool plt_dust_dens = param.getPlot(plIDnd) && (!data_pos_dd_list.empty()); //to
    // do if dust denity is possible
    bool plt_gas_temp = (data_pos_tg != MAX_UINT) && param.isInPlotList(GRIDgas_temp);
    //bool plt_dust_temp = (!data_pos_dt_list.empty()) && param.isInPlotList(GRIDdust_temp);
    bool plt_abs_ini = (!data_pos_abs_ini_list.empty()) && param.isInPlotList(GRIDabs_ini);
    bool plt_mag = (data_pos_mx != MAX_UINT) && param.isInPlotList(GRIDgas_dens);

    plt_mag = (data_pos_mx != MAX_UINT) && (data_pos_my != MAX_UINT) && (data_pos_my != MAX_UINT) &&
              param.isInPlotList(GRIDmx) && param.isInPlotList(GRIDmy) && param.isInPlotList(GRIDmz);

    plt_vel = (data_pos_vx != MAX_UINT) && (data_pos_vy != MAX_UINT) && (data_pos_vz != MAX_UINT) &&
              param.isInPlotList(GRIDvx) && param.isInPlotList(GRIDvy) && param.isInPlotList(GRIDvz);

    plt_delta = plt_gas_temp && plt_mag && (!data_pos_dt_list.empty());
    //plt_larm = plt_gas_temp && plt_mag && (!data_pos_dt_list.empty());
    plt_mach = plt_vel && plt_gas_temp;

    ullong per_counter = 0, per_max = bins * bins;
    string dens_filename = path + "gas_density.am";
    string dtemp_filename = path + "dust_temp.am";
    string gtemp_filename = path + "gas_temp.am";
    string magvec_filename = path + "mag_vec_field.am";
    string velvec_filename = path + "vel_vec_field.am";
    string magf_filename = path + "mag_field.am";
    string velf_filename = path + "vel_field.am";
    string a_filename = path + "aalig.am";
    string disr_filename = path + "adisr.am";
    string max_disr_filename = path + "max_adisr.am";
    string param_modif_filename = path + "param_modif.am";
    string barnet_low_lower_filename = path + "abar_low_lower.am";
    string barnet_low_upper_filename = path + "abar_low_upper.am"; 
    string barnet_high_lower_filename = path + "abar_high_lower.am";
    string barnet_high_upper_filename = path + "abar_high_upper.am";
    string adg_lower_filename = path + "adg_lower.am";
    string adg_upper_filename = path + "adg_upper.am";
    string adg_10_lower_filename = path + "adg_10_lower.am";
    string adg_10_upper_filename = path + "adg_10_upper.am";
    string abs_ini_filename = path + "abs_ini.am";
    string amaxJB_Lar_filename = path + "amaxJB_Lar.am";
    string akrat_lowJ_filename = path + "akrat_lowJ.am";
    string akrat_highJ_filename = path + "akrat_highJ.am";
    string d_filename = path + "delta.am";

    ofstream dens_writer, rat_writer, amaxJB_Lar_writer;
    ofstream akrat_lowJ_writer, akrat_highJ_writer;
    ofstream disr_writer, max_disr_writer, param_modif_writer;
    ofstream barnet_low_lower_writer, barnet_low_upper_writer;
    ofstream barnet_high_lower_writer, barnet_high_upper_writer;
    ofstream adg_lower_writer, adg_upper_writer, adg_10_lower_writer, adg_10_upper_writer;
    ofstream gas_writer, dust_writer;
    ofstream abs_ini_writer;
    ofstream magvec_writer, magf_writer;
    ofstream velvec_writer, velf_writer;
    ofstream d_writer;

    if(plt_gas_dens)
    {
        dens_writer.open(dens_filename.c_str(), ios::out);

        if(dens_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << dens_filename << endl;
            return false;
        }
    }

    if(plt_gas_temp)
    {
        gas_writer.open(gtemp_filename.c_str(), ios::out);

        if(gas_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << gtemp_filename << endl;
            return false;
        }
    }

    if(plt_dust_temp)
    {
        dust_writer.open(dtemp_filename.c_str(), ios::out);
        if(dust_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << dtemp_filename << endl;
            return false;
        }
    }

    if(plt_rat)
    {
        rat_writer.open(a_filename.c_str(), ios::out);

        if(rat_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << a_filename << endl;
            return false;
        }
    }

    if(plt_disr)
    {
        disr_writer.open(disr_filename.c_str(), ios::out);

        if(disr_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << disr_filename << endl;
            return false;
        }
    }

    if(plt_max_disr)
    {
        max_disr_writer.open(max_disr_filename.c_str(), ios::out);

        if(max_disr_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << max_disr_filename << endl;
            return false;
        }
    }

    if(plt_param_modif)
    {
        param_modif_writer.open(param_modif_filename.c_str(), ios::out);

        if(param_modif_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << param_modif_filename << endl;
            return false;
        }
    }
    
 
    if(plt_barnet_low_lower)
    {
        barnet_low_lower_writer.open(barnet_low_lower_filename.c_str(), ios::out);

        if(barnet_low_lower_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << barnet_low_lower_filename << endl;
            return false;
        }
    }

    if(plt_barnet_low_upper)
    {
        barnet_low_upper_writer.open(barnet_low_upper_filename.c_str(), ios::out);

        if(barnet_low_upper_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << barnet_low_upper_filename << endl;
            return false;
        }
    }
    
    
    if(plt_barnet_high_lower)
    {
        barnet_high_lower_writer.open(barnet_high_lower_filename.c_str(), ios::out);

        if(barnet_high_lower_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << barnet_high_lower_filename << endl;
            return false;
        }
    }

    if(plt_barnet_high_upper)
    {
        barnet_high_upper_writer.open(barnet_high_upper_filename.c_str(), ios::out);

        if(barnet_high_upper_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << barnet_high_upper_filename << endl;
            return false;
        }
    }
    
    if(plt_dg_lower)
    {
        adg_lower_writer.open(adg_lower_filename.c_str(), ios::out);

        if(adg_lower_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << adg_lower_filename << endl;
            return false;
        }
    }

    if(plt_dg_upper)
    {
        adg_upper_writer.open(adg_upper_filename.c_str(), ios::out);

        if(adg_upper_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << adg_upper_filename << endl;
            return false;
        }
    }
    
    if(plt_dg_10_lower)
    {
        adg_10_lower_writer.open(adg_10_lower_filename.c_str(), ios::out);

        if(adg_10_lower_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << adg_10_lower_filename << endl;
            return false;
        }
    }

    if(plt_dg_10_upper)
    {
        adg_10_upper_writer.open(adg_10_upper_filename.c_str(), ios::out);

        if(adg_10_upper_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << adg_10_upper_filename << endl;
            return false;
        }
    }
    
    if(plt_abs_ini)
    {
        abs_ini_writer.open(abs_ini_filename.c_str(), ios::out);
        if(abs_ini_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << abs_ini_filename << endl;
            return false;
        }
    }
    
    
    if(plt_amaxJB_Lar)
    {
        amaxJB_Lar_writer.open(amaxJB_Lar_filename.c_str(), ios::out);

        if(amaxJB_Lar_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << amaxJB_Lar_filename << endl;
            return false;
        }
    }

    if(plt_akrat_lowJ)
    {
        akrat_lowJ_writer.open(akrat_lowJ_filename.c_str(), ios::out);

        if(akrat_lowJ_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << akrat_lowJ_filename << endl;
            return false;
        }
    }

    if(plt_akrat_highJ)
    {
        akrat_highJ_writer.open(akrat_highJ_filename.c_str(), ios::out);

        if(akrat_lowJ_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << akrat_highJ_filename << endl;
            return false;
        }
    }
    
    if(plt_delta)
    {
        d_writer.open(d_filename.c_str(), ios::out);

        if(d_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << d_filename << endl;
            return false;
        }
    }

    if(plt_mag)
    {
        magvec_writer.open(magvec_filename.c_str(), ios::out);

        if(magvec_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << magvec_filename << endl;
            return false;
        }

        magf_writer.open(magf_filename.c_str(), ios::out);

        if(magf_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << magf_filename << endl;
            return false;
        }
    }

    if(plt_vel)
    {
        velvec_writer.open(velvec_filename.c_str(), ios::out);

        if(velvec_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << velvec_filename << endl;
            return false;
        }

        velf_writer.open(velf_filename.c_str(), ios::out);

        if(velf_writer.fail())
        {
            cout << "\nERROR: Cannot write to:\n " << velf_filename << endl;
            return false;
        }
    }

    stringstream point_header, vec_header;
    point_header.str("");
    vec_header.str("");

    int b_limit = int(bins) / 2;
    double xyz_step = max_len / double(bins);

    double off_xyz = 0.5 * xyz_step;
    photon_package * pp = new photon_package;

    point_header << "# AmiraMesh 3D ASCII 2.0\n" << endl;
    point_header << "define Lattice " << bins << " " << bins << " " << bins;
    point_header << "\tParameters {" << endl;
    point_header << "Content \"" << bins << "x" << bins << "x" << bins << " float, uniform coordinates\","
                 << endl;
    point_header << "BoundingBox ";
    point_header << 0 << " "; // float(cell_oc_root->getXmin())
    point_header << 1 << " "; // float(cell_oc_root->getXmax())
    point_header << 0 << " "; // float(cell_oc_root->getXmin())
    point_header << 1 << " "; // float(cell_oc_root->getXmax())
    point_header << 0 << " "; // float(cell_oc_root->getXmin())
    point_header << 1 << " "; // float(cell_oc_root->getXmax())
    point_header << "," << endl;
    point_header << " CoordType \"uniform\"" << endl;
    point_header << "}" << endl;
    point_header << "Lattice { float Data } @1" << endl;
    point_header << "# Data section follows" << endl;
    point_header << "@1" << endl;

    vec_header << "# AmiraMesh 3D ASCII 2.0\n" << endl;
    vec_header << "define Lattice " << bins << " " << bins << " " << bins;
    vec_header << "\tParameters {" << endl;
    vec_header << "Content \"" << bins << "x" << bins << "x" << bins << " float[3], uniform coordinates\","
               << endl;
    vec_header << "BoundingBox ";
    vec_header << 0 << " "; // float(cell_oc_root->getXmin())
    vec_header << 1 << " "; // float(cell_oc_root->getXmax())
    vec_header << 0 << " "; // float(cell_oc_root->getXmin())
    vec_header << 1 << " "; // float(cell_oc_root->getXmax())
    vec_header << 0 << " "; // float(cell_oc_root->getXmin())
    vec_header << 1 << " "; // float(cell_oc_root->getXmax())
    vec_header << "," << endl;
    vec_header << " CoordType \"uniform\"" << endl;
    vec_header << "}" << endl;
    vec_header << "Lattice { float[3] Data } @1" << endl;
    vec_header << "# Data section follows" << endl;
    vec_header << "@1" << endl;

    dens_writer << "# " << (min_len) << " " << (max_len);
    dens_writer << point_header.str();
    gas_writer << point_header.str();
    dust_writer << point_header.str();
    abs_ini_writer << point_header.str();
    magf_writer << point_header.str();
    velf_writer << point_header.str();

    magvec_writer << vec_header.str();
    velvec_writer << vec_header.str();

    rat_writer << point_header.str();
    disr_writer << point_header.str();
    max_disr_writer << point_header.str();
    param_modif_writer << point_header.str();
    barnet_low_lower_writer << point_header.str();
    barnet_low_upper_writer << point_header.str(); 
    barnet_high_lower_writer << point_header.str();
    barnet_high_upper_writer << point_header.str();
    adg_lower_writer << point_header.str();
    adg_upper_writer << point_header.str();
    adg_10_lower_writer << point_header.str();
    adg_10_upper_writer << point_header.str();
    akrat_lowJ_writer << point_header.str();
    akrat_highJ_writer << point_header.str();
    d_writer << point_header.str();
    amaxJB_Lar_writer << point_header.str();

    Vector3D mag_field, vel_field;

    for(int z = -b_limit; z <= b_limit; z++)
    {
        if(z == 0)
            continue;

        for(int y = -b_limit; y <= b_limit; y++)
        {
            if(y == 0)
                continue;

            for(int x = -b_limit; x <= b_limit; x++)
            {
                if(x == 0)
                    continue;

                double sgx = CMathFunctions::sgn(x);
                double sgy = CMathFunctions::sgn(y);
                double sgz = CMathFunctions::sgn(z);
                double tx = double(x) * xyz_step - sgx * off_xyz;
                double ty = double(y) * xyz_step - sgy * off_xyz;
                double tz = double(z) * xyz_step - sgz * off_xyz;

                pp->setPosition(Vector3D(tx, ty, tz));

                if(!positionPhotonInGrid(pp))
                {
                    dens_writer << float(log10(min_gas_dens)) << endl;
                    gas_writer << uint(0) << endl;
                    dust_writer << uint(0) << endl;
                    magf_writer << float(log10(min_mag)) << endl;
                    velf_writer << float(log10(min_vel)) << endl;

                    magvec_writer << uint(0) << " " << uint(0) << " " << uint(0) << endl;
                    velvec_writer << uint(0) << " " << uint(0) << " " << uint(0) << endl;

                    d_writer << uint(0) << endl;

                    rat_writer << uint(0) << endl;
                    disr_writer << uint(0) << endl;
                    max_disr_writer << uint(0) << endl;
                    param_modif_writer << uint(0) << endl;
                    barnet_low_lower_writer << uint(0) << endl;
                    barnet_low_upper_writer << uint(0) << endl; 
                    barnet_high_lower_writer << uint(0) << endl;
                    barnet_high_upper_writer << uint(0) << endl;
                    adg_lower_writer << uint(0) << endl;
                    adg_upper_writer << uint(0) << endl;
                    adg_10_lower_writer << uint(0) << endl;
                    adg_10_upper_writer << uint(0) << endl;
                    abs_ini_writer << uint(0) << endl;
                    amaxJB_Lar_writer << uint(0) << endl;
                    akrat_lowJ_writer << uint(0) << endl;
                    akrat_highJ_writer << uint(0) << endl;
                }
                else
                {
                    dens_writer << float(log10(getGasDensity(pp))) << endl;
                    gas_writer << float((getGasTemperature(pp))) << endl;
                    dust_writer << float((getDustTemperature(pp))) << endl;
                    abs_ini_writer << float((getQBOffset(pp))) << endl;

                    if(plt_mag)
                    {
                        mag_field = getMagField(pp);
                        magf_writer << float(log10(mag_field.length())) << endl;
                        mag_field.normalize();
                        // mag_field *= max_len;

                        magvec_writer << float(mag_field.X()) << " " << float(mag_field.Y()) << " "
                                      << float(mag_field.Z()) << endl;
                    }

                    if(plt_vel)
                    {
                        vel_field = getVelocityField(pp);
                        velf_writer << float(log10(vel_field.length())) << endl;

                        vel_field.normalize();
                        // vel_field *= max_len;

                        velvec_writer << float(vel_field.X()) << " " << float(vel_field.Y()) << " "
                                      << float(vel_field.Z()) << endl;
                    }

                    if(plt_rat)
                        rat_writer << float(log10(getAlignedRadius(pp, 0))) << endl;

                    if(plt_disr)
                        disr_writer << float(log10(getDisruptRadius(pp, 0))) << endl;

                    if(plt_max_disr)
                        max_disr_writer << float(log10(getMaxDisruptRadius(pp, 0))) << endl;

                    if(plt_param_modif)
                        param_modif_writer << float((getSizeParamModify(pp, 0))) << endl;
                     
                    if(plt_barnet_low_lower)
                        barnet_low_lower_writer << float((getBarnetLowLowerRadius(pp, 0))) << endl;
                        
                    if(plt_barnet_low_upper)
                        barnet_low_upper_writer << float((getBarnetLowUpperRadius(pp, 0))) << endl;                   
                        
                    if(plt_barnet_high_lower)
                        barnet_high_lower_writer << float((getBarnetHighLowerRadius(pp, 0))) << endl;
                        
                    if(plt_barnet_high_upper)
                        barnet_high_upper_writer << float((getBarnetHighUpperRadius(pp, 0))) << endl;
                        
                    if(plt_dg_lower)
                        adg_lower_writer << float((getDGLowerRadius(pp, 0))) << endl;
                        
                    if(plt_dg_upper)
                        adg_upper_writer << float((getDGUpperRadius(pp, 0))) << endl;
                        
                    if(plt_dg_10_lower)
                        adg_10_lower_writer << float((getDG10LowerRadius(pp, 0))) << endl;
                        
                    if(plt_dg_10_upper)
                        adg_10_upper_writer << float((getDG10UpperRadius(pp, 0))) << endl;

                    if(plt_amaxJB_Lar)
                        amaxJB_Lar_writer << float((getMaxAlignedRadius(pp, 0))) << endl;

                    if(plt_akrat_lowJ)
                        akrat_lowJ_writer << float((getkRATlowJRadius(pp, 0))) << endl;

                    if(plt_akrat_highJ)
                        akrat_highJ_writer << float((getkRAThighJRadius(pp, 0))) << endl;
                        
                    if(plt_delta)
                    {
                        double field = getMagField(pp).length();
                        double Td = getDustTemperature(pp);
                        double Tg = getGasTemperature(pp);
                        double dens = getGasDensity(pp);
                        double delta = CMathFunctions::calc_delta(field, Td, Tg, dens);
                        d_writer << float(log10(delta)) << endl;
                    }
                }
            }

            per_counter++;
            if(per_counter % 49 == 0)
                cout << " -> Writing AMIRA files:     " << 100.0 * float(per_counter) / float(per_max)
                     << " [%]                  \r";
        }
    }

    dens_writer.close();
    gas_writer.close();
    dust_writer.close();
    magvec_writer.close();
    magf_writer.close();

    velvec_writer.close();
    velf_writer.close();
    rat_writer.close();
    disr_writer.close();
    max_disr_writer.close();
    param_modif_writer.close();
    barnet_low_lower_writer.close();
    barnet_low_upper_writer.close(); 
    barnet_high_lower_writer.close();
    barnet_high_upper_writer.close();
    adg_lower_writer.close();
    adg_upper_writer.close();
    adg_10_lower_writer.close();
    adg_10_upper_writer.close();
    amaxJB_Lar_writer.close();
    abs_ini_writer.close();
    akrat_lowJ_writer.close();
    akrat_highJ_writer.close();
    d_writer.close();

    cout << "- Writing AMIRA files                  : done" << endl;

    delete pp;
    return true;
}

bool CGridBasic::writePhysicalParameter(string data_path)
{
    string x_filename = data_path + "physical_parameter_x.dat";
    string y_filename = data_path + "physical_parameter_y.dat";
    string z_filename = data_path + "physical_parameter_z.dat";

    ofstream x_writer, y_writer, z_writer;

    x_writer.open(x_filename.c_str(), ios::out);
    y_writer.open(y_filename.c_str(), ios::out);
    z_writer.open(z_filename.c_str(), ios::out);

    if(x_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n " << x_filename << endl;
        return false;
    }

    if(y_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n " << y_filename << endl;
        return false;
    }

    if(z_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n " << z_filename << endl;
        return false;
    }

    x_writer.precision(8);
    x_writer << scientific;

    y_writer.precision(8);
    y_writer << scientific;

    z_writer.precision(8);
    z_writer << scientific;

    cout << " -> Writing lines: 0.0 [%]                   \r" << flush;

    photon_package * pp = new photon_package;

    // along z
    pp->setPosition(Vector3D(0, 0, 2.0 * max_len));
    pp->setDirection(Vector3D(0.0001, 0.0001, -1.00001).normalized());
    findStartingPoint(pp);

    z_writer << "d\tng\tTg\tTd\tmx\tmy\tmz\tvx\tvy\tvz" << endl;


    while(next(pp))
    {
        double pos = pp->getPosition().Z();
        double dens = getGasDensity(pp);
        double Tg = 0;
        double Td = 0;
        double mx = 0, my = 0, mz = 0;
        double vx = 0, vy = 0, vz = 0;

        if(data_pos_tg != MAX_UINT)
            Tg = getGasTemperature(pp);

        if(!data_pos_dt_list.empty())
            Td = getDustTemperature(pp);

        if(data_pos_mx != MAX_UINT)
        {
            mx = getMagField(pp).X();
            my = getMagField(pp).Y();
            mz = getMagField(pp).Z();
        }

        if(data_pos_vx != MAX_UINT)
        {
            vx = getVelocityField(pp).X();
            vy = getVelocityField(pp).Y();
            vz = getVelocityField(pp).Z();
        }
            
        z_writer << pos << "\t" << dens << "\t" << Tg << "\t" << Td << "\t" << mx << "\t" << my << "\t" << mz
                 << "\t" << vx << "\t" << vy << "\t" << vz << endl;
 
    }

    cout << " -> Writing lines: 33.3 [%]                  \r" << flush;

    pp->setPosition(Vector3D(0, 2.0 * max_len, 0));
    pp->setDirection(Vector3D(0.0001, -1.00001, 0.0001).normalized());
    findStartingPoint(pp);

    y_writer << "d\tng\tTg\tTd\tmx\tmy\tmz\tvx\tvy\tvz" << endl;

    while(next(pp))
    {
        double pos = pp->getPosition().Y();
        double dens = getGasDensity(pp);
        double Tg = 0;
        double Td = 0;
        double mx = 0, my = 0, mz = 0;
        double vx = 0, vy = 0, vz = 0;

        if(data_pos_tg != MAX_UINT)
            Tg = getGasTemperature(pp);

        if(!data_pos_dt_list.empty())
            Td = getDustTemperature(pp);

        if(data_pos_mx != MAX_UINT)
        {
            mx = getMagField(pp).X();
            my = getMagField(pp).Y();
            mz = getMagField(pp).Z();
        }

        if(data_pos_vx != MAX_UINT)
        {
            vx = getVelocityField(pp).X();
            vy = getVelocityField(pp).Y();
            vz = getVelocityField(pp).Z();
        }
 
        y_writer << pos << "\t" << dens << "\t" << Tg << "\t" << Td << "\t" << mx << "\t" << my << "\t" << mz
                 << "\t" << vx << "\t" << vy << "\t" << vz << endl;
    }

    cout << " -> Writing lines: 66.6 [%]                   \r" << flush;

    pp->setPosition(Vector3D(2.0 * max_len, 0, 0));
    pp->setDirection(Vector3D(-1.00001, 0.0001, 0.0001).normalized());
    findStartingPoint(pp);

    x_writer << "d\tng\tTg\tTd\tmx\tmy\tmz\tvx\tvy\tvz" << endl;
    while(next(pp))
    {
        double pos = pp->getPosition().X();
        double dens = getGasDensity(pp);
        double Tg = 0;
        double Td = 0;
        double mx = 0, my = 0, mz = 0;
        double vx = 0, vy = 0, vz = 0;

        if(data_pos_tg != MAX_UINT)
            Tg = getGasTemperature(pp);

        if(!data_pos_dt_list.empty())
            Td = getDustTemperature(pp);

        if(data_pos_mx != MAX_UINT)
        {
            mx = getMagField(pp).X();
            my = getMagField(pp).Y();
            mz = getMagField(pp).Z();
        }

        if(data_pos_vx != MAX_UINT)
        {
            vx = getVelocityField(pp).X();
            vy = getVelocityField(pp).Y();
            vz = getVelocityField(pp).Z();
        }
 
        x_writer << pos << "\t" << dens << "\t" << Tg << "\t" << Td << "\t" << mx << "\t" << my << "\t" << mz
                 << "\t" << vx << "\t" << vy << "\t" << vz << endl;
    }

    x_writer.close();
    y_writer.close();
    z_writer.close();

    delete pp;

    cout << "- Writing lines                 : done" << endl;
    return true;
}



bool CGridBasic::writeSizeParameter(string data_path)
{
    string x_filename = data_path + "grain_size_x.dat";
    string y_filename = data_path + "grain_size_y.dat";
    string z_filename = data_path + "grain_size_z.dat";

    ofstream x_writer, y_writer, z_writer;

    x_writer.open(x_filename.c_str(), ios::out);
    y_writer.open(y_filename.c_str(), ios::out);
    z_writer.open(z_filename.c_str(), ios::out);

    if(x_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n " << x_filename << endl;
        return false;
    }

    if(y_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n " << y_filename << endl;
        return false;
    }

    if(z_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n " << z_filename << endl;
        return false;
    }

    x_writer.precision(8);
    x_writer << scientific;

    y_writer.precision(8);
    y_writer << scientific;

    z_writer.precision(8);
    z_writer << scientific;

    cout << " -> Writing lines: 0.0 [%]                   \r" << flush;

    photon_package * pp = new photon_package;

    // along z
    pp->setPosition(Vector3D(0, 0, 2.0 * max_len));
    pp->setDirection(Vector3D(0.0001, 0.0001, -1.00001).normalized());
    findStartingPoint(pp);

    z_writer << "d\ta_alg\tamaxJB_Lar\tadisr\tmax_adisr\ta_barnet_low_J_lower\ta_barnet_low_J_upper\ta_barnet_high_J_lower\ta_barnet_high_J_upper\ta_dg_lower\ta_dg_upper\ta_dg_10_lower\ta_dg_10_upper\takrat_lowJ\takrat_highJ"<< endl;


    while(next(pp))
    {
        double pos = pp->getPosition().Z();
        double a_alg = 0;
        double amaxJB_Lar = 0;
        double adisr = 0;
        double max_adisr = 0;
        double abar_low_lower = 0;
        double abar_low_upper = 0; 
        double abar_high_lower = 0;
        double abar_high_upper = 0;
        double adg_lower = 0;
        double adg_upper = 0;
        double adg_10_lower = 0;
        double adg_10_upper = 0;
        double akrat_lowJ = 0;
        double akrat_highJ = 0;

        if(!data_pos_aalg_list.empty())
            a_alg = getAlignedRadius(pp, 0);
        
        if(!data_pos_amaxJB_Lar_list.empty())
            amaxJB_Lar = getMaxAlignedRadius(pp, 0);

        if(!data_pos_adisr_list.empty())
            adisr = getDisruptRadius(pp, 0);

        if(!data_pos_max_adisr_list.empty())
            max_adisr = getMaxDisruptRadius(pp, 0);
            
        if(!data_pos_barnet_low_J_lower_list.empty())
            abar_low_lower = getBarnetLowLowerRadius(pp, 0);
            
        if(!data_pos_barnet_low_J_upper_list.empty())
            abar_low_upper = getBarnetLowUpperRadius(pp, 0); 
            
        if(!data_pos_barnet_high_J_lower_list.empty())
            abar_high_lower = getBarnetHighLowerRadius(pp, 0);
            
        if(!data_pos_barnet_high_J_upper_list.empty())
            abar_high_upper = getBarnetHighUpperRadius(pp, 0);
            
        if(!data_pos_dg_lower_list.empty())
            adg_lower = getDGLowerRadius(pp, 0);
            
        if(!data_pos_dg_upper_list.empty())
            adg_upper = getDGUpperRadius(pp, 0);
            
        if(!data_pos_dg_10_lower_list.empty())
            adg_10_lower = getDG10LowerRadius(pp, 0);
            
        if(!data_pos_dg_10_upper_list.empty())
            adg_10_upper = getDG10UpperRadius(pp, 0);

        if(!data_pos_akrat_lowJ_list.empty())
            akrat_lowJ = getkRATlowJRadius(pp, 0);

        if(!data_pos_akrat_highJ_list.empty())
            akrat_highJ = getkRAThighJRadius(pp, 0);
           
        z_writer << pos << "\t" << a_alg << "\t" << amaxJB_Lar << "\t" << adisr << "\t" << max_adisr << "\t" 
        << abar_low_lower << "\t" << abar_low_upper << "\t" << abar_high_lower << "\t" << abar_high_upper << 
        "\t" << adg_lower << "\t" << adg_upper << "\t" << adg_10_lower << "\t" << adg_10_upper  << "\t" << 
        akrat_lowJ << "\t" << akrat_highJ << endl;
 
    }

    cout << " -> Writing lines: 33.3 [%]                  \r" << flush;

    pp->setPosition(Vector3D(0, 2.0 * max_len, 0));
    pp->setDirection(Vector3D(0.0001, -1.00001, 0.0001).normalized());
    findStartingPoint(pp);
    
    
    y_writer << "d\ta_alg\tamaxJB_Lar\tadisr\tmax_adisr\ta_barnet_low_J_lower\ta_barnet_low_J_upper\ta_barnet_high_J_lower\ta_barnet_high_J_upper\ta_dg_lower\ta_dg_upper\ta_dg_10_lower\ta_dg_10_upper\rakrat_lowJ\rakrat_highJ" << endl;

    while(next(pp))
    {
        double pos = pp->getPosition().Y();
        double a_alg = 0;
        double amaxJB_Lar = 0;
        double adisr = 0;
        double max_adisr = 0;
        double abar_low_lower = 0;
        double abar_low_upper = 0; 
        double abar_high_lower = 0;
        double abar_high_upper = 0;
        double adg_lower = 0;
        double adg_upper = 0;
        double adg_10_lower = 0;
        double adg_10_upper = 0;
        double akrat_lowJ = 0;
        double akrat_highJ = 0;

        if(!data_pos_aalg_list.empty())
            a_alg = getAlignedRadius(pp, 0);

        if(!data_pos_amaxJB_Lar_list.empty())
            amaxJB_Lar = getMaxAlignedRadius(pp, 0);
            
        if(!data_pos_adisr_list.empty())
            adisr = getDisruptRadius(pp, 0);

        if(!data_pos_max_adisr_list.empty())
            max_adisr = getMaxDisruptRadius(pp, 0);
            
        if(!data_pos_barnet_low_J_lower_list.empty())
            abar_low_lower = getBarnetLowLowerRadius(pp, 0);
            
        if(!data_pos_barnet_low_J_upper_list.empty())
            abar_low_upper = getBarnetLowUpperRadius(pp, 0); 
            
        if(!data_pos_barnet_high_J_lower_list.empty())
            abar_high_lower = getBarnetHighLowerRadius(pp, 0);
            
        if(!data_pos_barnet_high_J_upper_list.empty())
            abar_high_upper = getBarnetHighUpperRadius(pp, 0);
            
        if(!data_pos_dg_lower_list.empty())
            adg_lower = getDGLowerRadius(pp, 0);
            
        if(!data_pos_dg_upper_list.empty())
            adg_upper = getDGUpperRadius(pp, 0);
            
        if(!data_pos_dg_10_lower_list.empty())
            adg_10_lower = getDG10LowerRadius(pp, 0);
            
        if(!data_pos_dg_10_upper_list.empty())
            adg_10_upper = getDG10UpperRadius(pp, 0);

        if(!data_pos_akrat_lowJ_list.empty())
            akrat_lowJ = getkRATlowJRadius(pp, 0);

        if(!data_pos_akrat_highJ_list.empty())
            akrat_highJ = getkRAThighJRadius(pp, 0);
 
        y_writer << pos << "\t" << a_alg << "\t" << amaxJB_Lar << "\t" << adisr << "\t" << max_adisr << "\t" 
        << abar_low_lower << "\t" << abar_low_upper << "\t" << abar_high_lower << "\t" << abar_high_upper << 
        "\t" << adg_lower << "\t" << adg_upper << "\t" << adg_10_lower << "\t" << adg_10_upper  << "\t" << 
        akrat_lowJ << "\t" << akrat_highJ << endl;
    }

    cout << " -> Writing lines: 66.6 [%]                   \r" << flush;

    pp->setPosition(Vector3D(2.0 * max_len, 0, 0));
    pp->setDirection(Vector3D(-1.00001, 0.0001, 0.0001).normalized());
    findStartingPoint(pp);

    x_writer << "d\ta_alg\tamaxJB_Lar\tadisr\tmax_adisr\ta_barnet_low_J_lower\ta_barnet_low_J_upper\ta_barnet_high_J_lower\ta_barnet_high_J_upper\ta_dg_lower\ta_dg_upper\ta_dg_10_lower\ta_dg_10_upper\takrat_lowJ\tlt_akrat_highJ" << endl;
    while(next(pp))
    {
        double pos = pp->getPosition().X();
        double a_alg = 0;
        double adisr = 0;
        double amaxJB_Lar = 0;
        double max_adisr = 0;
        double abar_low_lower = 0;
        double abar_low_upper = 0; 
        double abar_high_lower = 0;
        double abar_high_upper = 0;
        double adg_lower = 0;
        double adg_upper = 0;
        double adg_10_lower = 0;
        double adg_10_upper = 0;
        double akrat_lowJ = 0;
        double akrat_highJ = 0;

        if(!data_pos_aalg_list.empty())
            a_alg = getAlignedRadius(pp, 0);
            
        if(!data_pos_amaxJB_Lar_list.empty())
            amaxJB_Lar = getMaxAlignedRadius(pp, 0);

        if(!data_pos_adisr_list.empty())
            adisr = getDisruptRadius(pp, 0);

        if(!data_pos_max_adisr_list.empty())
            max_adisr = getMaxDisruptRadius(pp, 0);
            
        if(!data_pos_barnet_low_J_lower_list.empty())
            abar_low_lower = getBarnetLowLowerRadius(pp, 0);
            
        if(!data_pos_barnet_low_J_upper_list.empty())
            abar_low_upper = getBarnetLowUpperRadius(pp, 0); 
            
        if(!data_pos_barnet_high_J_lower_list.empty())
            abar_high_lower = getBarnetHighLowerRadius(pp, 0);
            
        if(!data_pos_barnet_high_J_upper_list.empty())
            abar_high_upper = getBarnetHighUpperRadius(pp, 0);
            
        if(!data_pos_dg_lower_list.empty())
            adg_lower = getDGLowerRadius(pp, 0);
            
        if(!data_pos_dg_upper_list.empty())
            adg_upper = getDGUpperRadius(pp, 0);
            
        if(!data_pos_dg_10_lower_list.empty())
            adg_10_lower = getDG10LowerRadius(pp, 0);
            
        if(!data_pos_dg_10_upper_list.empty())
            adg_10_upper = getDG10UpperRadius(pp, 0);

        if(!data_pos_akrat_lowJ_list.empty())
            akrat_lowJ = getkRATlowJRadius(pp, 0);

        if(!data_pos_akrat_highJ_list.empty())
            akrat_highJ = getkRAThighJRadius(pp, 0);
            
        x_writer << pos << "\t" << a_alg << "\t" << amaxJB_Lar << "\t" << adisr << "\t" << max_adisr << "\t" 
        << abar_low_lower << "\t" << abar_low_upper << "\t" << abar_high_lower << "\t" << abar_high_upper << 
        "\t" << adg_lower << "\t" << adg_upper << "\t" << adg_10_lower << "\t" << adg_10_upper  << "\t" << 
        akrat_lowJ << "\t" << akrat_highJ << endl; 
    }

    x_writer.close();
    y_writer.close();
    z_writer.close();

    delete pp;

    cout << "- Writing lines                 : done" << endl;
    return true;
}

bool CGridBasic::writeMidplaneFits(string data_path, parameters & param, uint bins, bool all)
{
    bool res = true;

    if(bins == 0)
        return res;

    int cmd = param.getCommand();

    if(all)
    {
        plt_gas_dens = (!data_pos_gd_list.empty()) && param.isInPlotList(GRIDgas_dens);
        plt_dust_dens = (!data_pos_dd_list.empty()) && param.isInPlotList(GRIDdust_dens);
        plt_gas_temp = (data_pos_tg != MAX_UINT) && param.isInPlotList(GRIDgas_temp);

        plt_mag = (data_pos_mx != MAX_UINT) && (data_pos_my != MAX_UINT) && (data_pos_my != MAX_UINT) &&
                  param.isInPlotList(GRIDmx) && param.isInPlotList(GRIDmy) && param.isInPlotList(GRIDmz);

        plt_vel = (data_pos_vx != MAX_UINT) && (data_pos_vy != MAX_UINT) && (data_pos_vz != MAX_UINT) &&
                  param.isInPlotList(GRIDvx) && param.isInPlotList(GRIDvy) && param.isInPlotList(GRIDvz);

        plt_delta = plt_gas_temp && plt_mag && (!data_pos_dt_list.empty());
        plt_mach = plt_vel && plt_gas_temp;

        plt_dust_id = (data_pos_id != MAX_UINT);
		plt_amin = (data_pos_amin != MAX_UINT) && param.isInPlotList(GRIDa_min);
        plt_amax = (data_pos_amax != MAX_UINT) && param.isInPlotList(GRIDa_max);
        plt_size_param = (data_pos_size_param != MAX_UINT) && param.isInPlotList(GRIDq);

        plt_n_th = (data_pos_n_th != MAX_UINT) && param.isInPlotList(GRIDn_th);
		plt_T_e = (data_pos_T_e != MAX_UINT) && param.isInPlotList(GRIDT_e);
        plt_n_cr = (data_pos_n_cr != MAX_UINT) && param.isInPlotList(GRIDn_cr);
        plt_g_min = (data_pos_g_min != MAX_UINT) && param.isInPlotList(GRIDg_min);
        plt_g_max = (data_pos_g_max != MAX_UINT) && param.isInPlotList(GRIDg_max);
        plt_p = (data_pos_p != MAX_UINT) && param.isInPlotList(GRIDp);
 
        if(cmd != CMD_RAT && cmd != CMD_TEMP_RAT && cmd != CMD_TEMP_RAT_DISR && cmd != CMD_RAT_DISR)
        {
            plt_rat = (!data_pos_aalg_list.empty()) && param.isInPlotList(GRIDa_alg);
            plt_amaxJB_Lar = (!data_pos_amaxJB_Lar_list.empty()) && param.isInPlotList(GRIDamaxJB_Lar);
            plt_avg_th = (data_pos_avg_th != MAX_UINT) && param.isInPlotList(GRIDavg_th);
            plt_avg_dir = (data_pos_avg_dir != MAX_UINT) && param.isInPlotList(GRIDavg_dir);
            
            if (param.getAligMRAT())
            {
    		    plt_barnet_low_lower = (!data_pos_barnet_low_J_lower_list.empty()) && param.isInPlotList(GRIDabar_low_lower);
    		    plt_barnet_low_upper = (!data_pos_barnet_low_J_upper_list.empty()) && param.isInPlotList(GRIDabar_low_upper); 
    		    plt_barnet_high_lower = (!data_pos_barnet_high_J_lower_list.empty()) && param.isInPlotList(GRIDabar_high_lower);
    		    plt_barnet_high_upper = (!data_pos_barnet_high_J_upper_list.empty()) && param.isInPlotList(GRIDabar_high_upper);
    		    plt_dg_lower = (!data_pos_dg_lower_list.empty()) && param.isInPlotList(GRIDadg_lower);
    		    plt_dg_upper = (!data_pos_dg_upper_list.empty()) && param.isInPlotList(GRIDadg_upper);
    		    plt_dg_10_lower = (!data_pos_dg_10_lower_list.empty()) && param.isInPlotList(GRIDadg_10_lower);
    		    plt_dg_10_upper = (!data_pos_dg_10_upper_list.empty()) && param.isInPlotList(GRIDadg_10_upper);
            }
            if (param.getAligkRAT())
            {
            	plt_akrat_lowJ = (!data_pos_akrat_lowJ_list.empty()) && param.isInPlotList(GRIDakrat_lowJ);
            	plt_akrat_highJ = (!data_pos_akrat_highJ_list.empty()) && param.isInPlotList(GRIDakrat_highJ);
            }
        }

        if(cmd != CMD_DISR && cmd != CMD_TEMP_DISR && cmd != CMD_TEMP_RAT_DISR && CMD_RAT_DISR)
        {
            plt_disr = (!data_pos_adisr_list.empty()) && param.isInPlotList(GRIDadisr);
            plt_max_disr = (!data_pos_max_adisr_list.empty()) && param.isInPlotList(GRIDadisr_max);
            plt_param_modif = (!data_pos_param_modif_list.empty()) && param.isInPlotList(GRIDparam_modif);
            plt_avg_th = (data_pos_avg_th != MAX_UINT) && param.isInPlotList(GRIDavg_th);
            plt_avg_dir = (data_pos_avg_dir != MAX_UINT) && param.isInPlotList(GRIDavg_dir);
        }

        //if(cmd != CMD_TEMP_DISR || cmd != CMD_DISR)
        //{
        //    plt_disr = (!data_pos_adisr_list.empty()) && param.isInPlotList(GRIDadisr);
        //    plt_max_disr = (!data_pos_max_adisr_list.empty()) && param.isInPlotList(GRIDadisr_max);
        //    plt_param_modif = (!data_pos_param_modif_list.empty()) && param.isInPlotList(GRIDparam_modif);
        //    plt_avg_th = (data_pos_avg_th != MAX_UINT) && param.isInPlotList(GRIDavg_th);
        //    plt_avg_dir = (data_pos_avg_dir != MAX_UINT) && param.isInPlotList(GRIDavg_dir);
        //}

        if(cmd != CMD_TEMP && cmd != CMD_TEMP_RAT && CMD_TEMP_DISR && CMD_TEMP_RAT_DISR)
	    {
            plt_dust_temp = (!data_pos_dt_list.empty()) && param.isInPlotList(GRIDdust_temp);
            //plt_abs_ini = (!data_pos_abs_ini_list.empty()) && param.isInPlotList(GRIDabs_ini);
 	    }
        if(getRadiationFieldAvailable())
        {
            switch(param.getWriteRadiationField())
            {
                default:
                    plt_u_rad = false;
                    plt_rad_field1 = false;
                    break;

                case 1:
                    plt_u_rad = (cmd == CMD_RAT || cmd == CMD_TEMP_RAT || cmd == CMD_DISR || cmd == CMD_TEMP_DISR || cmd == CMD_RAT_DISR || cmd == CMD_TEMP_RAT_DISR);
                    plt_rad_field1 = false;
                    break;

                case 2:
                    plt_u_rad = false;
                    plt_rad_field1 = true;
                    break;

                case 3:
                    plt_u_rad = false;
                    plt_rad_field1 = true;
                    nr_rad_field_comp = 4;
                    break;
            }

            if(param.getWriteGZero())
                plt_g_zero1 = true;
        }
    }
    else
    {
        plt_gas_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_mag = false;
        plt_vel = false;
        plt_rat = false;
        plt_delta = false;
        //plt_larm = false;
        plt_mach = false;
        plt_dust_id = false;
        plt_amin = false;
        plt_amax = false;
        plt_size_param = false;
        plt_rad_field1 = false;
        plt_g_zero1 = false;
        plt_u_rad = false;
        plt_n_th = false;
        plt_T_e = false;
        plt_n_cr = false;
        plt_g_min = false;
        plt_g_max = false;
        plt_p = false;

        plt_avg_th = false;
        plt_avg_dir = false;

        plt_disr = false;
        plt_max_disr = false;
        plt_param_modif = false;
        plt_barnet_low_lower = false;
        plt_barnet_low_upper = false; 
        plt_barnet_high_lower = false;
        plt_barnet_high_upper = false;
        plt_dg_lower = false;
        plt_dg_upper = false;
        plt_dg_10_lower = false;
        plt_dg_10_upper = false;
        //plt_abs_ini = false;
        plt_amaxJB_Lar = false;
        plt_akrat_lowJ = false;
        plt_akrat_highJ = false;


        if(cmd == CMD_TEMP || cmd == CMD_TEMP_RAT || cmd == CMD_TEMP_DISR || cmd == CMD_TEMP_RAT_DISR)
        {
            if(param.getAdjTgas() > 0)
                plt_gas_temp = param.isInPlotList(GRIDgas_temp);

            plt_dust_temp = param.isInPlotList(GRIDdust_temp);
            //plt_abs_ini = param.isInPlotList(GRIDabs_ini);
        }
        if(cmd == CMD_RAT || cmd == CMD_TEMP_RAT || cmd == CMD_RAT_DISR || cmd == CMD_TEMP_RAT_DISR)
        {
            plt_rat = param.isInPlotList(GRIDa_alg);
            plt_amaxJB_Lar  = param.isInPlotList(GRIDamaxJB_Lar);
            plt_avg_th = param.isInPlotList(GRIDavg_th);
            plt_avg_dir = param.isInPlotList(GRIDavg_dir);
            if (param.getAligMRAT())
            {
		    plt_barnet_low_lower  = param.isInPlotList(GRIDabar_low_lower);
		    plt_barnet_low_upper  = param.isInPlotList(GRIDabar_low_upper); 
		    plt_barnet_high_lower  = param.isInPlotList(GRIDabar_high_lower);
		    plt_barnet_high_upper  = param.isInPlotList(GRIDabar_high_upper);
		    plt_dg_lower  = param.isInPlotList(GRIDadg_lower);
		    plt_dg_upper  = param.isInPlotList(GRIDadg_upper);
		    plt_dg_10_lower  = param.isInPlotList(GRIDadg_10_lower);
		    plt_dg_10_upper  = param.isInPlotList(GRIDadg_10_upper);
 	    }
 	    if (param.getAligkRAT())
 	    {
            	plt_akrat_lowJ  = param.isInPlotList(GRIDakrat_lowJ);
            	plt_akrat_highJ  = param.isInPlotList(GRIDakrat_highJ);
            }
        }

        if(cmd == CMD_DISR || cmd == CMD_TEMP_DISR || cmd == CMD_RAT_DISR || cmd == CMD_TEMP_RAT_DISR)
        {
            plt_disr = param.isInPlotList(GRIDadisr);
            plt_max_disr = param.isInPlotList(GRIDadisr_max);
            plt_param_modif = param.isInPlotList(GRIDparam_modif);
            plt_avg_th = param.isInPlotList(GRIDavg_th);
            plt_avg_dir = param.isInPlotList(GRIDavg_dir);
        }

        //if(cmd == CMD_TEMP_DISR || cmd == CMD_DISR)
        //{
        //    plt_disr = param.isInPlotList(GRIDadisr);
        //    plt_max_disr = param.isInPlotList(GRIDadisr_max);
        //    plt_param_modif = param.isInPlotList(GRIDparam_modif);
        //    plt_avg_th = param.isInPlotList(GRIDavg_th);
        //    plt_avg_dir = param.isInPlotList(GRIDavg_dir);
        //}

        switch(param.getWriteRadiationField())
        {
            default:
                plt_u_rad = false;
                plt_rad_field1 = false;
                break;

            case 1:
                plt_u_rad = (cmd == CMD_RAT || cmd == CMD_TEMP_RAT || cmd == CMD_DISR || cmd == CMD_TEMP_DISR || cmd == CMD_RAT_DISR ||  cmd == CMD_TEMP_RAT_DISR);
                plt_rad_field1 = false;
                break;

            case 2:
                plt_u_rad = false;
                plt_rad_field1 = true;

                if(!spec_length_as_vector)
                    cout << "HINT: The full radiation field can only be saved if it was used by the "
                            "simulation\n"
                            "      (when saving the radiation field in the grid or calculating RATs)!"
                         << endl;
                break;

            case 3:
                plt_u_rad = false;
                plt_rad_field1 = true;
                nr_rad_field_comp = 4;
                break;
        }

        if(param.getWriteGZero())
            plt_g_zero1 = true;
    }

    uint nr_parameters = uint(plt_gas_dens) + uint(plt_dust_dens) + uint(plt_gas_temp) + uint(plt_dust_temp) +
                         4 * uint(plt_mag) + 4 * uint(plt_vel) + uint(plt_rat) + + uint(plt_amaxJB_Lar) +
                         uint(plt_delta) + uint(plt_mach) + uint(plt_dust_id) + uint(plt_rad_field1) * 
                         nr_rad_field_comp * WL_STEPS + uint(plt_g_zero1) + uint(plt_u_rad) + uint(plt_n_th) + 
                         uint(plt_T_e) + uint(plt_n_cr) + uint(plt_g_min) + uint(plt_g_max) + uint(plt_p) + 
                         uint(plt_avg_th) + uint(plt_avg_dir) + uint(plt_disr) + uint(plt_max_disr) + 
                         uint(plt_param_modif) + uint(plt_barnet_low_lower) + uint(plt_barnet_low_upper) + 
                         uint(plt_barnet_high_lower) + uint(plt_barnet_high_upper) + uint(plt_dg_lower) + 
                         uint(plt_dg_upper) + uint(plt_dg_10_lower) + uint(plt_dg_10_upper) + 
                         uint(plt_akrat_lowJ) + uint(plt_akrat_highJ);

    if(nr_parameters == 0)
        return res;

    if(plt_gas_dens)
        if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
            nr_parameters += nr_densities;
    if(plt_dust_dens)
        if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
            nr_parameters += nr_densities;
    if(plt_dust_temp)
        if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
            nr_parameters += nr_densities;
    //if(plt_abs_ini)
    //    if(nr_densities > 1 && data_pos_abs_ini_list.size() >= nr_densities)
    //        nr_parameters += nr_densities;
 

    long naxis = 4;
    long naxes[4] = { uint(bins), uint(bins), 3, nr_parameters };
    uint per_max = 3 * bins * bins;

    double max_midplane_len = (max_len / param.getMidplaneZoom());

    dlist midplane_3d_param = param.getMidplane3dParams();
    double z_step, off_z, shift_z = 0;
    uint plane_3d = 0;
    if(midplane_3d_param.size() == 4)
    {
        plane_3d = midplane_3d_param[0];

        if(midplane_3d_param[1] != 0)
        {
            naxes[2] = uint(midplane_3d_param[1]);
            per_max = bins * bins * midplane_3d_param[1];
        }
        else
        {
            naxes[2] = uint(bins);
            per_max = bins * bins * bins;
        }

        if(midplane_3d_param[2] != 0 || midplane_3d_param[3] != 0)
        {
            z_step = (midplane_3d_param[3] - midplane_3d_param[2]) / double(naxes[2]);
            off_z = 0.5 * z_step;
            shift_z = (midplane_3d_param[3] + midplane_3d_param[2]) / 2.0;
        }
        else
        {
            z_step = max_midplane_len / double(naxes[2]);
            off_z = 0.5 * z_step;
        }
    }
    else
    {
        z_step = max_midplane_len / double(bins);
        off_z = 0.5 * z_step;
    }

    double xy_step = max_midplane_len / double(bins);
    double off_xy = 0.5 * xy_step;
    int b_limit_z, b_limit_xy;

    if(naxes[2] % 2)
    {
        b_limit_z = (naxes[2] - 1) / 2;
        off_z = 0;
    }
    else
        b_limit_z = naxes[2] / 2;

    if(naxes[0] % 2)
    {
        b_limit_xy = (naxes[0] - 1) / 2;
        off_xy = 0;
    }
    else
        b_limit_xy = naxes[0] / 2;

    ullong per_counter = 0;

    auto_ptr<CCfits::FITS> pFits(0);
    // unique_ptr<CCfits::FITS> pFits;

    try
    {
        string path_out = data_path + "midplane.fits";
        if(midplane_3d_param.size() == 4)
            path_out = data_path + "midplane_3d.fits";
        remove(path_out.c_str());
        pFits.reset(new CCfits::FITS(path_out, DOUBLE_IMG, naxis, naxes));
    }
    catch(CCfits::FITS::CantCreate)
    {
        return false;
    }

    long nelements = bins * bins;

    valarray<double> array_gas_dens(nelements);
    valarray<double> array_dust_dens(nelements);
    valarray<double> array_gas_temp(nelements);
    valarray<double> array_dust_temp(nelements);
    valarray<double> array_rat(nelements);
    valarray<double> array_delta(nelements);
    valarray<double> array_mag(nelements);
    valarray<double> array_mag_x(nelements);
    valarray<double> array_mag_y(nelements);
    valarray<double> array_mag_z(nelements);
    valarray<double> array_vel(nelements);
    valarray<double> array_vel_x(nelements);
    valarray<double> array_vel_y(nelements);
    valarray<double> array_vel_z(nelements);
    //valarray<double> array_larm(nelements);
    valarray<double> array_mach(nelements);
    valarray<double> array_dust_mixture(nelements);
    valarray<double> array_amin(nelements);
    valarray<double> array_amax(nelements);
    valarray<double> array_size_param(nelements);
    valarray<double> array_rad_field(nelements);
    valarray<double> array_g_zero1(nelements);
    valarray<double> array_u_rad(nelements);
    valarray<double> array_n_th(nelements);
    valarray<double> array_T_e(nelements);
    valarray<double> array_n_cr(nelements);
    valarray<double> array_g_min(nelements);
    valarray<double> array_g_max(nelements);
    valarray<double> array_p(nelements);

    valarray<double> array_avg_th(nelements);
    valarray<double> array_avg_dir(nelements);
   

    valarray<double> array_disr(nelements);
    valarray<double> array_max_disr(nelements);
    valarray<double> array_param_modif(nelements);
    valarray<double> array_barnet_low_lower(nelements);
    valarray<double> array_barnet_low_upper(nelements); 
    valarray<double> array_barnet_high_lower(nelements);
    valarray<double> array_barnet_high_upper(nelements);
    valarray<double> array_dg_lower(nelements);
    valarray<double> array_dg_upper(nelements);
    valarray<double> array_dg_10_lower(nelements);
    valarray<double> array_dg_10_upper(nelements);
    //valarray<double> array_abs_ini(nelements);
    valarray<double> array_amaxJB_Lar(nelements);
    valarray<double> array_akrat_lowJ(nelements);
    valarray<double> array_akrat_highJ(nelements);
    
    if(plt_gas_dens)
    {
        buffer_gas_dens = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are
            // in the grid
            if(nr_densities > 1 && data_pos_gd_list.size() == nr_densities)
                buffer_gas_dens[i_cell] = new double[nr_densities + 1];
            else
                buffer_gas_dens[i_cell] = new double[nr_densities];
        }
    }
    
    if(plt_dust_dens)
    {
        buffer_dust_dens = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are
            // in the grid
            if(nr_densities > 1 && data_pos_dd_list.size() == nr_densities)
                buffer_dust_dens[i_cell] = new double[nr_densities + 1];
            else
                buffer_dust_dens[i_cell] = new double[nr_densities];
        }
    }
    
    if(plt_gas_temp)
        buffer_gas_temp = new double[nelements];
        
    if(plt_dust_temp)
    {
        buffer_dust_temp = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            // +1 for the average/sum of the quantity, but only if multiple quantities are
            // in the grid
            if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                buffer_dust_temp[i_cell] = new double[nr_densities + 1];
            else
                buffer_dust_temp[i_cell] = new double[nr_densities];
        }
    }
    
    if(plt_rat)
    {
        buffer_rat = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_rat[i_cell] = new double[data_pos_aalg_list.size()];
    }
    
    if(plt_delta)
        buffer_delta = new double[nelements];
        
    if(plt_mag)
    {
        buffer_mag = new double[nelements];
        buffer_mag_x = new double[nelements];
        buffer_mag_y = new double[nelements];
        buffer_mag_z = new double[nelements];
    }
    
    if(plt_vel)
    {
        buffer_vel = new double[nelements];
        buffer_vel_x = new double[nelements];
        buffer_vel_y = new double[nelements];
        buffer_vel_z = new double[nelements];
    }
 
    if(plt_mach)
        buffer_mach = new double[nelements];
        
    if(plt_dust_id)
        buffer_dust_mixture = new double[nelements];
        
    if(plt_amin)
        buffer_dust_amin = new double[nelements];
        
    if(plt_amax)
        buffer_dust_amax = new double[nelements];
        
    if(plt_size_param)
        buffer_dust_size_param = new double[nelements];
        
    if(plt_rad_field1)
    {
        buffer_rad_field = new double **[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            buffer_rad_field[i_cell] = new double *[WL_STEPS];
            for(int wID = 0; wID < WL_STEPS; wID++)
                buffer_rad_field[i_cell][wID] = new double[nr_rad_field_comp];
        }
    }
    
    if(plt_g_zero1)
        buffer_g_zero1 = new double[nelements];
        
    if(plt_u_rad)
        buffer_u_rad = new double[nelements];
        
    if(plt_n_th)
        buffer_n_th = new double[nelements];
        
    if(plt_T_e)
        buffer_T_e = new double[nelements];
        
    if(plt_n_cr)
        buffer_n_cr = new double[nelements];
        
    if(plt_g_min)
        buffer_g_min = new double[nelements];
        
    if(plt_g_max)
        buffer_g_max = new double[nelements];
        
    if(plt_p)
        buffer_p = new double[nelements];
        
    if(plt_avg_th)
        buffer_avg_th = new double[nelements];
        
    if(plt_avg_dir)
        buffer_avg_dir = new double[nelements];
        
      
    if(plt_disr)
    {
        buffer_disr = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_disr[i_cell] = new double[data_pos_adisr_list.size()];
    }
    
    if(plt_max_disr)
    {
        buffer_max_disr = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_max_disr[i_cell] = new double[data_pos_max_adisr_list.size()];
    }
    
    if(plt_param_modif)
    {
        buffer_param_modif = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_param_modif[i_cell] = new double[data_pos_param_modif_list.size()];
    }
    
    
    if(plt_barnet_low_lower)
    {
        buffer_barnet_low_lower = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_barnet_low_lower[i_cell] = new double[data_pos_barnet_low_J_lower_list.size()];
    }
    if(plt_barnet_low_upper)
    {
        buffer_barnet_low_upper = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_barnet_low_upper[i_cell] = new double[data_pos_barnet_low_J_upper_list.size()];
    } 
    
    
    
    if(plt_barnet_high_lower)
    {
        buffer_barnet_high_lower = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_barnet_high_lower[i_cell] = new double[data_pos_barnet_high_J_lower_list.size()];
    }
    if(plt_barnet_high_upper)
    {
        buffer_barnet_high_upper = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_barnet_high_upper[i_cell] = new double[data_pos_barnet_high_J_upper_list.size()];
    }
    
    
    
    if(plt_dg_lower)
    {
        buffer_dg_lower = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_dg_lower[i_cell] = new double[data_pos_dg_lower_list.size()];
    }
    if(plt_dg_upper)
    {
        buffer_dg_upper = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_dg_upper[i_cell] = new double[data_pos_dg_upper_list.size()];
    }
    
    
    if(plt_dg_10_lower)
    {
        buffer_dg_10_lower = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_dg_10_lower[i_cell] = new double[data_pos_dg_10_lower_list.size()];
    }
    if(plt_dg_10_upper)
    {
        buffer_dg_10_upper = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_dg_10_upper[i_cell] = new double[data_pos_dg_10_upper_list.size()];
    }
    
  
    //if(plt_abs_ini)
    //{
    //    buffer_abs_ini = new double *[nelements];
    //    for(int i_cell = 0; i_cell < nelements; i_cell++)
    //    {
    //        // +1 for the average/sum of the quantity, but only if multiple quantities are
    //        // in the grid
    //        if(nr_densities > 1 && data_pos_abs_ini_list.size() >= nr_densities)
    //            buffer_abs_ini[i_cell] = new double[nr_densities + 1];
    //        else
    //            buffer_abs_ini[i_cell] = new double[nr_densities];
    //    }
    //}
    
    if(plt_amaxJB_Lar)
    {
        buffer_amaxJB_Lar = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_amaxJB_Lar[i_cell] = new double[data_pos_amaxJB_Lar_list.size()];
    }

    if(plt_akrat_lowJ)
    {
        buffer_akrat_lowJ = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_akrat_lowJ[i_cell] = new double[data_pos_akrat_lowJ_list.size()];
    }

    if(plt_akrat_highJ)
    {
        buffer_akrat_highJ = new double *[nelements];
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            buffer_akrat_highJ[i_cell] = new double[data_pos_akrat_highJ_list.size()];
    }


    

    vector<long> fpixel(4);

    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = 0;

    if(midplane_3d_param.size() == 4)
    {
        for(int l = -b_limit_z; l <= b_limit_z; l++)
        {
            if(l == 0 && naxes[2] % 2 == 0)
                continue;

            fpixel[2]++;

#pragma omp parallel for schedule(dynamic)
            for(long i_cell = 0; i_cell < nelements; i_cell++)
            {
                int j = (i_cell % bins);

                int k = i_cell / bins - b_limit_xy;

                j -= b_limit_xy;

                if(bins % 2 == 0)
                    if(k > -1)
                        k++;

                if(bins % 2 == 0)
                    if(j > -1)
                        j++;

                double tx, ty, tz;
                setPlaneParameter(plane_3d, xy_step, off_xy, z_step, off_z, shift_z, j, k, l, tx, ty, tz);

                fillMidplaneBuffer(tx, ty, tz, i_cell, param);

                per_counter++;
                if(per_counter % 220 == 0)
                {
#pragma omp critical
                    {
                        cout << " -> Writing 3D midplane file: "
                             << 100.0 * float(per_counter) / float(per_max) << " [%]                 \r";
                    }
                }
            }

            fpixel[3] = 0;
            
            if(plt_gas_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_dens[i_cell] = buffer_gas_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_dens);

                if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_gas_dens[i_cell] = buffer_gas_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_gas_dens);
                    }
            }
            
            if(plt_dust_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_dens[i_cell] = buffer_dust_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_dens);

                if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_dens[i_cell] = buffer_dust_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_dens);
                    }
            }
            
            if(plt_gas_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_temp[i_cell] = buffer_gas_temp[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_temp);
            }
            
            if(plt_dust_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_temp[i_cell] = buffer_dust_temp[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_temp);

                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_temp[i_cell] = buffer_dust_temp[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_temp);
                    }
            }
            
            if(plt_rat)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_rat[i_cell] = buffer_rat[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_rat);
                }
            }
            
            if(plt_delta)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_delta[i_cell] = buffer_delta[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_delta);
            }
            
            if(plt_mag)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_mag[i_cell] = buffer_mag[i_cell];
                    array_mag_x[i_cell] = buffer_mag_x[i_cell];
                    array_mag_y[i_cell] = buffer_mag_y[i_cell];
                    array_mag_z[i_cell] = buffer_mag_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_z);
            }
            
            if(plt_vel)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_vel[i_cell] = buffer_vel[i_cell];
                    array_vel_x[i_cell] = buffer_vel_x[i_cell];
                    array_vel_y[i_cell] = buffer_vel_y[i_cell];
                    array_vel_z[i_cell] = buffer_vel_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_z);
            }
 
            if(plt_mach)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_mach[i_cell] = buffer_mach[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mach);
            }
            
            if(plt_dust_id)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_mixture[i_cell] = buffer_dust_mixture[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_mixture);
            }
            
            if(plt_amin)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amin[i_cell] = buffer_dust_amin[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amin);
            }
            
            if(plt_amax)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amax[i_cell] = buffer_dust_amax[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amax);
            }
            
            if(plt_size_param)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_size_param[i_cell] = buffer_dust_size_param[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_size_param);
            }
            
            if(plt_rad_field1)
            {
                for(int i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                    for(int wID = 0; wID < WL_STEPS; wID++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_rad_field[i_cell] = buffer_rad_field[i_cell][wID][i_comp];

                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_rad_field);
                    }
            }
            
            if(plt_g_zero1)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_zero1[i_cell] = buffer_g_zero1[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_zero1);
            }
            
            if(plt_u_rad)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_u_rad[i_cell] = buffer_u_rad[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_u_rad);
            }
            
            if(plt_n_th)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_th[i_cell] = buffer_n_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_th);
            }
            
            if(plt_T_e)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_T_e[i_cell] = buffer_T_e[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_T_e);
            }
            
            if(plt_n_cr)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_cr[i_cell] = buffer_n_cr[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_cr);
            }
            
            if(plt_g_min)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_min[i_cell] = buffer_g_min[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_min);
            }
            
            if(buffer_g_max)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_max[i_cell] = buffer_g_max[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_max);
            }
            
            if(plt_p)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_p[i_cell] = buffer_p[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_p);
            }
            
            if(plt_avg_th)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_th[i_cell] = buffer_avg_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_th);
            }
            
            if(plt_avg_dir)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_dir[i_cell] = buffer_avg_dir[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_dir);
            }

            if(plt_disr)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_disr[i_cell] = buffer_disr[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_disr);
                }
            }

            if(plt_max_disr)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_max_disr[i_cell] = buffer_max_disr[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_max_disr);
                }
            }

            if(plt_param_modif)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_param_modif[i_cell] = buffer_param_modif[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_param_modif);
                }
            }
            
            if(plt_barnet_low_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_low_lower[i_cell] = buffer_barnet_low_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_low_lower);
                }
            }
            
            if(plt_barnet_low_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_low_upper[i_cell] = buffer_barnet_low_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_low_upper);
                }
            } 
            
            if(plt_barnet_high_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_high_lower[i_cell] = buffer_barnet_high_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_high_lower);
                }
            }
            
            if(plt_barnet_high_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_high_upper[i_cell] = buffer_barnet_high_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_high_upper);
                }
            }
            
            
            if(plt_dg_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_lower[i_cell] = buffer_dg_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_lower);
                }
            }
            
            if(plt_dg_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_upper[i_cell] = buffer_dg_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_upper);
                }
            }
            
            
            if(plt_dg_10_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_10_lower[i_cell] = buffer_dg_10_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_10_lower);
                }
            }
            
            if(plt_dg_10_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_10_upper[i_cell] = buffer_dg_10_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_10_upper);
                }
            }
            
            
            //if(plt_abs_ini)
            //{
            //    for(int i_cell = 0; i_cell < nelements; i_cell++)
            //        array_abs_ini[i_cell] = buffer_abs_ini[i_cell][0];
            //    fpixel[3]++;
            //    pFits->pHDU().write(fpixel, nelements, array_dust_temp);

            //    if(nr_densities > 1 && data_pos_abs_ini_list.size() >= nr_densities)
            //    {
            //        for(uint i_density = 0; i_density < nr_densities; i_density++)
            //        {
            //            for(int i_cell = 0; i_cell < nelements; i_cell++)
            //                array_abs_ini[i_cell] = buffer_abs_ini[i_cell][i_density + 1];
            //            fpixel[3]++;
            //            pFits->pHDU().write(fpixel, nelements, array_dust_temp);
            //        }
            // 	}
            //}
            
            if(plt_amaxJB_Lar)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_amaxJB_Lar[i_cell] = buffer_amaxJB_Lar[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_amaxJB_Lar);
                }
            }

            if(plt_akrat_lowJ)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_akrat_lowJ[i_cell] = buffer_akrat_lowJ[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_akrat_lowJ);
                }
            }

            if(plt_akrat_highJ)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_akrat_highJ[i_cell] = buffer_akrat_highJ[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_akrat_highJ);
                }
            }
            
        }
    }
    else
    {
        for(int i = 1; i <= 3; i++)
        {
            fpixel[2] = i;

#pragma omp parallel for schedule(dynamic)
            for(long i_cell = 0; i_cell < nelements; i_cell++)
            {
                int j = (i_cell % bins);

                int k = i_cell / bins - b_limit_xy;

                j -= b_limit_xy;

                if(bins % 2 == 0)
                    if(k > -1)
                        k++;

                if(bins % 2 == 0)
                    if(j > -1)
                        j++;

                double tx, ty, tz;

                setPlaneParameter(i, xy_step, off_xy, 0, 0, 0, j, k, 0, tx, ty, tz);

                fillMidplaneBuffer(tx, ty, tz, i_cell, param);

                per_counter++;
                if(per_counter % 220 == 0)
                {
#pragma omp critical
                    {
                        cout << " -> Writing midplane files: " << 100.0 * float(per_counter) / float(per_max)
                             << " [%]             \r";
                    }
                }
            }

            fpixel[3] = 0;
            
            if(plt_gas_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_dens[i_cell] = buffer_gas_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_dens);

                if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_gas_dens[i_cell] = buffer_gas_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_gas_dens);
                    }
            }
            
            if(plt_dust_dens)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_dens[i_cell] = buffer_dust_dens[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_dens);

                if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_dens[i_cell] = buffer_dust_dens[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_dens);
                    }
            }
            
            if(plt_gas_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_gas_temp[i_cell] = buffer_gas_temp[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_gas_temp);
            }
            
            if(plt_dust_temp)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_temp[i_cell] = buffer_dust_temp[i_cell][0];
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_temp);

                if(nr_densities > 1 && data_pos_dt_list.size() >= nr_densities)
                    for(uint i_density = 0; i_density < nr_densities; i_density++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_dust_temp[i_cell] = buffer_dust_temp[i_cell][i_density + 1];
                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_dust_temp);
                    }
            }
            
            if(plt_rat)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_rat[i_cell] = buffer_rat[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_rat);
                }
            }
            
            if(plt_delta)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_delta[i_cell] = buffer_delta[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_delta);
            }
            
            if(plt_mag)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_mag[i_cell] = buffer_mag[i_cell];
                    array_mag_x[i_cell] = buffer_mag_x[i_cell];
                    array_mag_y[i_cell] = buffer_mag_y[i_cell];
                    array_mag_z[i_cell] = buffer_mag_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mag_z);
            }
            
            if(plt_vel)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                {
                    array_vel[i_cell] = buffer_vel[i_cell];
                    array_vel_x[i_cell] = buffer_vel_x[i_cell];
                    array_vel_y[i_cell] = buffer_vel_y[i_cell];
                    array_vel_z[i_cell] = buffer_vel_z[i_cell];
                }

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_x);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_y);
                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_vel_z);
            }

            if(plt_mach)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_mach[i_cell] = buffer_mach[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_mach);
            }
            
            if(plt_dust_id)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_dust_mixture[i_cell] = buffer_dust_mixture[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_dust_mixture);
            }
            
            if(plt_amin)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amin[i_cell] = buffer_dust_amin[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amin);
            }
            
            if(plt_amax)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_amax[i_cell] = buffer_dust_amax[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_amax);
            }
            
            if(plt_size_param)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_size_param[i_cell] = buffer_dust_size_param[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_size_param);
            }
            
            if(plt_rad_field1)
            {
                for(int i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
                    for(int wID = 0; wID < WL_STEPS; wID++)
                    {
                        for(int i_cell = 0; i_cell < nelements; i_cell++)
                            array_rad_field[i_cell] = buffer_rad_field[i_cell][wID][i_comp];

                        fpixel[3]++;
                        pFits->pHDU().write(fpixel, nelements, array_rad_field);
                    }
            }
            
            if(plt_g_zero1)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_zero1[i_cell] = buffer_g_zero1[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_zero1);
            }
            
            if(plt_u_rad)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_u_rad[i_cell] = buffer_u_rad[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_u_rad);
            }
            
            if(plt_n_th)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_th[i_cell] = buffer_n_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_th);
            }
            
            if(plt_T_e)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_T_e[i_cell] = buffer_T_e[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_T_e);
            }
            
            if(plt_n_cr)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_n_cr[i_cell] = buffer_n_cr[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_n_cr);
            }
            
            if(plt_g_min)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_min[i_cell] = buffer_g_min[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_min);
            }
            
            if(buffer_g_max)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_g_max[i_cell] = buffer_g_max[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_g_max);
            }
            
            if(plt_p)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_p[i_cell] = buffer_p[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_p);
            }
            
            if(plt_avg_th)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_th[i_cell] = buffer_avg_th[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_th);
            }
            
            if(plt_avg_dir)
            {
                for(int i_cell = 0; i_cell < nelements; i_cell++)
                    array_avg_dir[i_cell] = buffer_avg_dir[i_cell];

                fpixel[3]++;
                pFits->pHDU().write(fpixel, nelements, array_avg_dir);
            }

            if(plt_disr)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_disr[i_cell] = buffer_disr[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_disr);
                }
            }


            if(plt_max_disr)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_max_disr[i_cell] = buffer_max_disr[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_max_disr);
                }
            }

            if(plt_param_modif)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_param_modif[i_cell] = buffer_param_modif[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_param_modif);
                }
            }
            
 
            if(plt_barnet_low_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_low_lower[i_cell] = buffer_barnet_low_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_low_lower);
                }
            }
            
            if(plt_barnet_low_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_low_upper[i_cell] = buffer_barnet_low_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_low_upper);
                }
            }
            
                        
            if(plt_barnet_high_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_high_lower[i_cell] = buffer_barnet_high_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_high_lower);
                }
            }
            
            if(plt_barnet_high_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_barnet_high_upper[i_cell] = buffer_barnet_high_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_barnet_high_upper);
                }
            }


            if(plt_dg_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_lower[i_cell] = buffer_dg_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_lower);
                }
            }
            
            if(plt_dg_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_upper[i_cell] = buffer_dg_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_upper);
                }
            }
            
            if(plt_dg_10_lower)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_10_lower[i_cell] = buffer_dg_10_lower[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_10_lower);
                }
            }
            
            if(plt_dg_10_upper)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_dg_10_upper[i_cell] = buffer_dg_10_upper[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_dg_10_upper);
                }
            }
            
                        
            //if(plt_abs_ini)
            //{
            //    for(int i_cell = 0; i_cell < nelements; i_cell++)
            //        array_abs_ini[i_cell] = buffer_abs_ini[i_cell][0];
            //    fpixel[3]++;
            //    pFits->pHDU().write(fpixel, nelements, array_abs_ini);

            //    if(nr_densities > 1 && data_pos_abs_ini_list.size() >= nr_densities)
            //    {
            //        for(uint i_density = 0; i_density < nr_densities; i_density++)
            //        {
            //            for(int i_cell = 0; i_cell < nelements; i_cell++)
            //                array_abs_ini[i_cell] = buffer_abs_ini[i_cell][i_density + 1];
            //            fpixel[3]++;
            //            pFits->pHDU().write(fpixel, nelements, array_abs_ini);
            //        }
            //    }
            //}
            
            if(plt_amaxJB_Lar)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_amaxJB_Lar[i_cell] = buffer_amaxJB_Lar[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_amaxJB_Lar);
                }
            }

            if(plt_akrat_lowJ)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_akrat_lowJ[i_cell] = buffer_akrat_lowJ[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_akrat_lowJ);
                }
            }

            if(plt_akrat_highJ)
            {
                for(uint i_density = 0; i_density < nr_densities; i_density++)
                {
                    for(int i_cell = 0; i_cell < nelements; i_cell++)
                        array_akrat_highJ[i_cell] = buffer_akrat_highJ[i_cell][i_density];

                    fpixel[3]++;
                    pFits->pHDU().write(fpixel, nelements, array_akrat_highJ);
                }
            }
        }
    }

    double bin_width = max_midplane_len / bins;
    double first_pix_val = -max_midplane_len / 2.0 + (bin_width / 2.0);

    // Grid
    pFits->pHDU().addKey("CTYPE1", "PARAM", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1", first_pix_val, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1", bin_width, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1", "m", "unit of axis 1");

    // Alternatively as AU grid
    pFits->pHDU().addKey("CTYPE1B", "PARAM", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1B", first_pix_val / con_AU, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1B", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1B", bin_width / con_AU, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1B", "AU", "unit of axis 1");

    // Alternatively as pc grid
    pFits->pHDU().addKey("CTYPE1C", "PARAM", "type of unit 1");
    pFits->pHDU().addKey("CRVAL1C", first_pix_val / con_pc, "value of axis 1");
    pFits->pHDU().addKey("CRPIX1C", 1, "pixel where CRVAL1 is defined ");
    pFits->pHDU().addKey("CDELT1C", bin_width / con_pc, "delta of axis 1");
    pFits->pHDU().addKey("CUNIT1C", "pc", "unit of axis 1");

    // Grid
    pFits->pHDU().addKey("CTYPE2", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2", first_pix_val, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2", bin_width, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2", "m", "unit of axis 2");

    // Alternatively as AU grid
    pFits->pHDU().addKey("CTYPE2B", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2B", first_pix_val / con_AU, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2B", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2B", bin_width / con_AU, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2B", "AU", "unit of axis 2");

    // Alternatively as pc grid
    pFits->pHDU().addKey("CTYPE2C", "PARAM", "type of unit 2");
    pFits->pHDU().addKey("CRVAL2C", first_pix_val / con_pc, "value of axis 2");
    pFits->pHDU().addKey("CRPIX2C", 1, "pixel where CRVAL2 is defined ");
    pFits->pHDU().addKey("CDELT2C", bin_width / con_pc, "delta of axis 2");
    pFits->pHDU().addKey("CUNIT2C", "pc", "unit of axis 2");
    if(midplane_3d_param.size() == 4)
    {
        double bin_width_z = z_step;
        double first_pix_val_z = shift_z - (z_step * double(naxes[2])) / 2.0 + (bin_width_z / 2.0);

        // Grid
        pFits->pHDU().addKey("CTYPE3", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3", first_pix_val_z, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3", bin_width_z, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3", "m", "unit of axis 3");

        // Alternatively as AU grid
        pFits->pHDU().addKey("CTYPE3B", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3B", first_pix_val_z / con_AU, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3B", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3B", bin_width_z / con_AU, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3B", "AU", "unit of axis 3");

        // Alternatively as pc grid
        pFits->pHDU().addKey("CTYPE3C", "PARAM", "type of unit 3");
        pFits->pHDU().addKey("CRVAL3C", first_pix_val_z / con_pc, "value of axis 3");
        pFits->pHDU().addKey("CRPIX3C", 1, "pixel where CRVAL3 is defined ");
        pFits->pHDU().addKey("CDELT3C", bin_width_z / con_pc, "delta of axis 3");
        pFits->pHDU().addKey("CUNIT3C", "pc", "unit of axis 3");

        // Quantities
        pFits->pHDU().addKey("CTYPE4", "PARAM", "type of unit 4");
        pFits->pHDU().addKey("CRVAL4", 1, "value of axis 4");
        pFits->pHDU().addKey("CRPIX4", 1, "pixel where CRVAL4 is defined ");
        pFits->pHDU().addKey("CDELT4", 1, "delta of axis 4");
        pFits->pHDU().addKey("CUNIT4", "see MIDPLANEX", "unit of axis 4");
    }

    uint counter = 0;
    char str_1[1024];
    char str_2[1024];
    if(plt_gas_dens)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        if(nr_densities > 1 && data_pos_gd_list.size() >= nr_densities)
        {
            if(gas_is_mass_density)
                pFits->pHDU().addKey(str_1, "total_gas_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "total_gas_number_density [m^-3]", str_2);
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3;
                if(gas_is_mass_density)
                    str_3 = getDensityString("gas_mass_density_%i [kg/m^3]", i_density);
                else
                    str_3 = getDensityString("gas_number_density_%i [m^-3]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        }
        else
        {
            if(gas_is_mass_density)
                pFits->pHDU().addKey(str_1, "gas_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "gas_number_density [m^-3]", str_2);
        }
    }
    if(plt_dust_dens)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        if(nr_densities > 1 && data_pos_dd_list.size() >= nr_densities)
        {
            if(dust_is_mass_density)
                pFits->pHDU().addKey(str_1, "total_dust_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "total_dust_number_density [m^-3]", str_2);
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3;
                if(dust_is_mass_density)
                    str_3 = getDensityString("dust_mass_density_%i [kg/m^3]", i_density);
                else
                    str_3 = getDensityString("dust_number_density_%i [m^-3]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        }
        else
        {
            if(dust_is_mass_density)
                pFits->pHDU().addKey(str_1, "dust_mass_density [kg/m^3]", str_2);
            else
                pFits->pHDU().addKey(str_1, "dust_number_density [m^-3]", str_2);
        }
    }
    if(plt_gas_temp)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "gas_temperature [K]", str_2);
    }
    if(plt_dust_temp)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        if(nr_densities > 1)
        {
            pFits->pHDU().addKey(str_1, "average_dust_temperature [K]", str_2);
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("dust_temperature_%i [K]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        }
        else
            pFits->pHDU().addKey(str_1, "dust_temperature [K]", str_2);
    }
    if(plt_rat)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("rat_aalig_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "rat_aalig [m]", str_2);
        }
    }
    if(plt_delta)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "delta [m]", str_2);
    }
    if(plt_mag)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_total [T]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_x [T]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_y [T]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mag_z [T]", str_2);
    }
    if(plt_vel)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_total [m/s]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_x [m/s]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_y [m/s]", str_2);
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "vel_z [m/s]", str_2);
    }
    if(plt_mach)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "mach number", str_2);
    }
    if(plt_dust_id)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "dust_mixture [index]", str_2);
    }
    if(plt_amin)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "a_min [m]", str_2);
    }
    if(plt_amax)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "a_max [m]", str_2);
    }
    if(plt_size_param)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "size param [value]", str_2);
    }
    if(plt_rad_field1)
    {
        for(int i_comp = 0; i_comp < nr_rad_field_comp; i_comp++)
        {
            for(int wID = 0; wID < WL_STEPS; wID++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                char str_3[1024];
                switch(i_comp)
                {
                    default:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;

                    case 1:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field_x [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field_x [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;

                    case 2:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field_y [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field_y [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;

                    case 3:
#ifdef WINDOWS
                        sprintf_s(str_3, "rad_field_z [W/m/m^2] (%.3e [m])", wl_list[wID]);
#else
                        sprintf(str_3, "rad_field_z [W/m/m^2] (%.3e [m])", wl_list[wID]);
#endif
                        break;
                }
                pFits->pHDU().addKey(str_1, string(str_3), str_2);
            }
        }
    }
    if(plt_g_zero1)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "G_0 (dustem)", str_2);
    }
    if(plt_u_rad)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "u_rad/u_isrf", str_2);
    }
    if(plt_n_th)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "therm_el_density [m^-3]", str_2);
    }
    if(plt_T_e)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "electron temperature [K]", str_2);
    }
    if(plt_n_cr)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "cr_el_density [m^-3]", str_2);
    }
    if(plt_g_min)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "gamma_min", str_2);
    }
    if(plt_g_max)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "gamma_max", str_2);
    }
    if(plt_p)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "syn_p", str_2);
    }
    if(plt_avg_th)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "avg. RAT cos(theta)", str_2);
    }
    if(plt_avg_dir)
    {
        counter++;
        updateMidplaneString(str_1, str_2, counter);
        pFits->pHDU().addKey(str_1, "avg. RAT aniso. (gamma)", str_2);
    }

    if(plt_disr)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("adisr_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "adisr [m]", str_2);
        }
    }

    if(plt_max_disr)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("max_adisr_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "max_adisr [m]", str_2);
        }
    }

    if(plt_param_modif)
    {
        if(nr_densities > 1)
        {
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("param_modify_%i", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "param_modify", str_2);
        }
    }
    
    if(plt_barnet_low_lower)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("amin_aJ_lowJ_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "amin_aJ_lowJ [m]", str_2);
        }
    }
    
    if(plt_barnet_low_upper)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("amax_aJ_lowJ_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "amax_aJ_lowJ [m]", str_2);
        }
    }
 

    if(plt_barnet_high_lower)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("amin_aJ_highJ_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "amin_aJ_highJ [m]", str_2);
        }
    }
    
    if(plt_barnet_high_upper)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("amax_aJ_highJ_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "amax_aJ_highJ [m]", str_2);
        }
    }
    
    
    
    if(plt_dg_lower)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("aminJB_DG_0.5%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "aminJB_DG_0.5 [m]", str_2);
        }
    }
    
    if(plt_dg_upper)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("amaxJB_DG_0.5_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "amaxJB_DG_0.5 [m]", str_2);
        }
    }
    
    if(plt_dg_10_lower)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("aminJB_DG_1_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "aminJB_DG_1[m]", str_2);
        }
    }
    
    if(plt_dg_10_upper)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("amaxJB_DG_1_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "amaxJB_DG_1 [m]", str_2);
        }
    }
    
    //if(plt_abs_ini)
    //{
    //    counter++;
    //    updateMidplaneString(str_1, str_2, counter);
    //    if(nr_densities > 1)
    //    {
    //        pFits->pHDU().addKey(str_1, "initial_absorption_rate", str_2);
    //        for(uint i_density = 1; i_density <= nr_densities; i_density++)
    //        {
    //            counter++;
    //            updateMidplaneString(str_1, str_2, counter);
    //            string str_3 = getDensityString("initial_absorption_rate_%i", i_density);
    //            pFits->pHDU().addKey(str_1, str_3, str_2);
    //        }
    //    }
    //    else
    //        pFits->pHDU().addKey(str_1, "initial_absorption_rate", str_2);
    //}
    
    if(plt_amaxJB_Lar)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("amaxJB_Lar_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "amaxJB_Lar [m]", str_2);
        }
    }

    if(plt_akrat_lowJ)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("akrat_lowJ_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "akrat_lowJ [m]", str_2);
        }
    }

    if(plt_akrat_highJ)
    {
        if(nr_densities > 1)
            for(uint i_density = 1; i_density <= nr_densities; i_density++)
            {
                counter++;
                updateMidplaneString(str_1, str_2, counter);
                string str_3 = getDensityString("akrat_highJ_%i [m]", i_density);
                pFits->pHDU().addKey(str_1, str_3, str_2);
            }
        else
        {
            counter++;
            updateMidplaneString(str_1, str_2, counter);
            pFits->pHDU().addKey(str_1, "akrat_highJ [m]", str_2);
        }
    }

    cout << CLR_LINE;
    cout << "Memory cleanup of the plotting arrays ...     \r" << flush;
    // Free memory of pointer arrays
    if(plt_gas_dens)
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_gas_dens[i_cell];
    if(plt_dust_dens)
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_dust_dens[i_cell];
    if(plt_gas_temp)
        delete[] buffer_gas_temp;
    if(plt_dust_temp)
        for(int i_cell = 0; i_cell < nelements; i_cell++)
            delete[] buffer_dust_temp[i_cell];
    if(plt_rat)
        delete[] buffer_rat;
    if(plt_delta)
        delete[] buffer_delta;
    if(plt_mag)
    {
        delete[] buffer_mag;
        delete[] buffer_mag_x;
        delete[] buffer_mag_y;
        delete[] buffer_mag_z;
    }
    if(plt_vel)
    {
        delete[] buffer_vel;
        delete[] buffer_vel_x;
        delete[] buffer_vel_y;
        delete[] buffer_vel_z;
    }
    //if(plt_larm)
    //    delete[] buffer_larm;
    if(plt_mach)
        delete[] buffer_mach;
    if(plt_dust_id)
        delete[] buffer_dust_mixture;
    if(plt_amin)
        delete[] buffer_dust_amin;
    if(plt_amax)
        delete[] buffer_dust_amax;
    if(plt_size_param)
        delete[] buffer_dust_size_param;
    if(plt_rad_field1)
    {
        for(int i_cell = 0; i_cell < nelements; i_cell++)
        {
            for(uint wID = 0; wID < WL_STEPS; wID++)
                delete[] buffer_rad_field[i_cell][wID];
            delete[] buffer_rad_field[i_cell];
        }
        delete[] buffer_rad_field;
    }
    if(plt_g_zero1)
        delete[] buffer_g_zero1;
    if(plt_u_rad)
        delete[] buffer_u_rad;
    if(plt_n_th)
        delete[] buffer_n_th;
    if(plt_T_e)
        delete[] buffer_T_e;
    if(plt_n_cr)
        delete[] buffer_n_cr;
    if(plt_g_min)
        delete[] buffer_g_min;
    if(plt_g_max)
        delete[] buffer_g_max;
    if(plt_p)
        delete[] buffer_p;
    if(plt_avg_th)
        delete[] buffer_avg_th;
    if(plt_avg_dir)
        delete[] buffer_avg_dir;

    if(plt_disr)
        delete[] buffer_disr;
    if(plt_max_disr)
        delete[] buffer_max_disr;
    if(plt_param_modif)
        delete[] buffer_param_modif;
    if(plt_barnet_low_lower)
        delete[] buffer_barnet_low_lower;
    if(plt_barnet_low_upper)
        delete[] buffer_barnet_low_upper;
    if(plt_barnet_high_lower)
        delete[] buffer_barnet_high_lower;
    if(plt_barnet_high_upper)
        delete[] buffer_barnet_high_upper;
    if(plt_dg_lower)
        delete[] buffer_dg_lower;
    if(plt_dg_upper)
        delete[] buffer_dg_upper;
    if(plt_dg_10_lower)
        delete[] buffer_dg_10_lower;
    if(plt_dg_10_upper)
        delete[] buffer_dg_10_upper;
    //if(plt_abs_ini)
    //{
    //    for(long i_cell = 0; i_cell < nelements; i_cell++)
    //        delete[] buffer_abs_ini[i_cell];
    //    delete[] buffer_abs_ini;
    //}
    if(plt_amaxJB_Lar)
        delete[] buffer_amaxJB_Lar;
    if(plt_akrat_lowJ)
        delete[] buffer_akrat_lowJ;
    if(plt_akrat_highJ)
        delete[] buffer_akrat_highJ;
        
    cout << CLR_LINE;
    cout << "- Writing of midplane files     : done" << endl;

    return res;
}

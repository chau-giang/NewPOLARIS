#include "Spherical.h"
#include "CommandParser.h"
#include "MathFunctions.h"
#include "typedefs.h"
#include <limits>

bool CGridSpherical::loadGridFromBinrayFile(parameters & param, uint _data_len)
{
    ushort tmpID, tmpOffset;
    string filename = param.getPathGrid();

    uint r_counter = 0;
    uint ph_counter = 0;
    uint th_counter = 0;

    line_counter = 0;
    char_counter = 0;

    turbulent_velocity = param.getTurbulentVelocity();

    ifstream bin_reader(filename.c_str(), ios::in | ios::binary);

    if(bin_reader.fail())
    {
        cout << "\nERROR: Cannot write to:\n Cannot load binary spherical grid file:" << endl;
        cout << filename << "\n" << endl;
        return false;
    }

    resetGridValues();

    max_cells = 0;

    line_counter = 1;
    char_counter = 0;
    float last_percentage = 0;

    bin_reader.read((char *)&tmpID, 2); // grid_ID
    bin_reader.read((char *)&tmpOffset, 2); // number of parameter in grid

    dataID = tmpID; // grid_ID
    data_offset = (uint)tmpOffset; // number of parameter in grid
    data_len = _data_len + data_offset; // Nparameter + parameter for use_energy_density?

    if(dataID == GRID_ID_SPH)
    {
        data_ids.resize(data_offset);

        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = 0;
            bin_reader.read((char *)&tmp_ids, 2);
            data_ids[i] = tmp_ids;
        }

        if(!setDataPositionsVariable())
            return false;
    }
    else
    {
        cout << "\nERROR: Cannot write to:\n A spherical grid requires an ID of \"" << GRID_ID_SPH << "\"!"
             << endl;
        return false;
    }

    uint tmp_data_offset = validateDataPositions(param);
    if(tmp_data_offset == MAX_UINT)
        return false;

    bin_reader.read((char *)&Rmin, 8);  
    bin_reader.read((char *)&Rmax, 8);
    bin_reader.read((char *)&N_r, 2);
    bin_reader.read((char *)&N_ph, 2);
    bin_reader.read((char *)&N_th, 2);
    bin_reader.read((char *)&log_factorR, 8);
    bin_reader.read((char *)&log_factorPh, 8);
    bin_reader.read((char *)&log_factorTh, 8);

    // Convert borders with conversion factors
    min_len = Rmin;
    Rmin *= conv_length_in_SI;
    Rmax *= conv_length_in_SI;

    total_volume = PIx4 * Rmax * Rmax * Rmax / 3.0; // volume of grid

    listR = new double[N_r + 1]; // array of cell point in radial direction
    listPh = new double[N_ph + 1]; // array of cell point in polar direction
    listTh = new double[N_th + 1]; // array of cell point in the azimuthal direction

    // --------------------------------------
    // ---------- Radial-direction ----------
    // --------------------------------------

    // Init radial cell border
    listR = new double[N_r + 1];
    if(log_factorR == 0)
    {
        // Allow user defined radius list, if log_factorR is zero

        // The global borders are already in the grid
        listR[0] = Rmin;
        listR[N_r] = Rmax;

        // Set the cell borders
        for(uint i_r = 1; i_r < N_r; i_r++)
        {
            // Read radial cell border position
            bin_reader.read((char *)&listR[i_r], 8);

            // Update radial position with conversion factors
            listR[i_r] *= conv_length_in_SI;
        }
    }
    else if(log_factorR == 1.0)
    {
        // Sinus shaped list, which emphasizes the middle rings
        CMathFunctions::SinList(Rmin, Rmax, listR, N_r + 1, log_factorR);
    }
    else if(log_factorR > 1.0)
    {
        // Exponentially increasing width of the cells in radial direction
        CMathFunctions::ExpList(Rmin, Rmax, listR, N_r + 1, log_factorR);
    }
    else
    {
        // Linear width of the cells in radial direction
        CMathFunctions::LinearList(Rmin, Rmax, listR, N_r + 1);
    }

    // -----------------------------------
    // ---------- Phi-direction ----------
    // -----------------------------------

    // Init phi cell border
    listPh = new double[N_ph + 1];
    if(log_factorPh == 0)
    {
        // Allow user defined phi list, if log_factorPh is zero

        // The global borders are already in the grid
        listPh[0] = 0;
        listPh[N_ph] = PIx2;

        // Set the cell borders
        for(uint i_ph = 1; i_ph < N_ph; i_ph++)
            bin_reader.read((char *)&listPh[i_ph], 8);
    }
    else
    {
        // Linear width of the cells in phi direction
        CMathFunctions::LinearList(0, PIx2, listPh, N_ph + 1);
    }

    // -------------------------------------
    // ---------- Theta-direction ----------
    // -------------------------------------
    if(log_factorTh == 0)
    {
        // The global borders are already in the grid
        listTh[0] = 0;
        listTh[N_th] = PI;

        // Read cell border in theta direction
        for(uint i_th = 1; i_th < N_th; i_th++)
            bin_reader.read((char *)&listTh[i_th], 8);
    }
    else if(log_factorTh == 1.0)
    {
        // Sinus shaped list, which emphasizes the midplane
        CMathFunctions::SinList(0, PI, listTh, N_th + 1, log_factorTh);
    }
    else if(log_factorTh > 1.0)
    {
        // Exponentially increasing width of the cells in z-direction (symmetrically)
        CMathFunctions::ExpListSym(0, PI, listTh, N_th + 1, log_factorTh);
    }
    else
    {
        // Linear width of the cells in theta direction
        CMathFunctions::LinearList(0, PI, listTh, N_th + 1);
    }

    // -----------------------------------------
    // ---------- Check of the limits ----------
    // -----------------------------------------
    if(Rmin <= 0)
    {
        cout << "\nERROR: Cannot write to:\n Inner radius (Rmin = " << Rmin << ") must be larger than zero!"
             << endl;
        return false;
    }

    if(Rmax <= 0)
    {
        cout << "\nERROR: Cannot write to:\n Outer radius (Rmax = " << Rmax << ") must be larger than zero!"
             << endl;
        return false;
    }

    if(Rmax <= Rmin)
    {
        cout << "\nERROR: Cannot write to:\n Outer radius (Rmax = " << Rmax
             << ") must be larger than inner radius (Rmin = " << Rmin << ")!" << endl;
        return false;
    }

    // Init grid cells
    grid_cells = new cell_sp ***[N_r];

    for(uint i_r = 0; i_r < N_r; i_r++) // with each value of r_i
    {
        grid_cells[i_r] = new cell_sp **[N_ph];

        cout << "Allocating memory for spherical grid cells: " << float(100.0 * double(i_r) / double(N_r))
             << "      \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++) // with each value of ph_i
        {
            grid_cells[i_r][i_ph] = new cell_sp *[N_th];

            for(uint i_th = 0; i_th < N_th; i_th++) // with each value of th_i
            {
                grid_cells[i_r][i_ph][i_th] = 0;
            }
        }
    }

    // Clear user output
    cout << CLR_LINE;

    max_cells = N_r * N_ph * N_th + 1;
    line_counter = -1;

    while(!bin_reader.eof())
    {
        line_counter++;

        if(r_counter < N_r)
        {
            double dr = listR[r_counter + 1] - listR[r_counter];

            if(dr == 0)
            {
                cout << "\nERROR: Cannot write to:\n No step size in r-direction of "
                        "spherical grid!"
                     << endl;
                return false;
            }

            if(dr < min_len)
                min_len = dr;
        }

        if(r_counter == 0)
        {
            double d;

            if(ph_counter < N_ph)
            {
                double dph = listPh[ph_counter + 1] - listPh[ph_counter];

                if(dph == 0)
                {
                    cout << "\nERROR: Cannot write to:\n No step size in phi-direction "
                            "of spherical grid!"
                         << endl;
                    return false;
                }

                if(N_ph > 2)
                {
                    d = 2 * Rmin * tan(dph / 2);
                    if(d < min_len)
                        min_len = d;
                }
            }

            if(th_counter < N_th)
            {
                double dth = listTh[th_counter + 1] - listTh[th_counter];

                if(dth == 0)
                {
                    cout << "\nERROR: Cannot write to:\n No step size in theta-direction "
                            "of spherical grid!"
                         << endl
                         << "\nHINT: Update of POLARIS v4.02 includes variable phi "
                            "spacing."
                         << endl
                         << "      Please look in the manual or use \"polaris-gen ... "
                            "--update\""
                         << endl;
                    return false;
                }

                d = 2 * Rmin * tan(dth / 2);
                if(d < min_len)
                    min_len = d;
            }
        }

        // Calculate percentage of total progress per source
        float percentage = 100.0 * double(line_counter) / double(max_cells);

        // Show only new percentage number if it changed
        if((percentage - last_percentage) > PERCENTAGE_STEP)
        {
            char_counter++;
            cout << "-> Loading spherical grid file: " << percentage << " [%]      \r" << flush;
            last_percentage = percentage;
        }

        if(th_counter == N_th)
        {
            ph_counter++;
            th_counter = 0;
        }

        if(ph_counter == N_ph)
        {
            r_counter++;
            ph_counter = 0;
        }

        if(r_counter > N_r)
            break;

        cell_sp * tmp_cell = 0;

        if(r_counter == N_r)
        {
            center_cell = new cell_sp;
            center_cell->setRID(MAX_UINT);
            center_cell->setPhID(MAX_UINT);
            center_cell->setThID(MAX_UINT);
            center_cell->resize(data_len + tmp_data_offset);
            tmp_cell = center_cell;
            tmp_cell->setID(0);
            r_counter++;
        }
        else
        {
            grid_cells[r_counter][ph_counter][th_counter] = new cell_sp;
            grid_cells[r_counter][ph_counter][th_counter]->setRID(r_counter);
            grid_cells[r_counter][ph_counter][th_counter]->setPhID(ph_counter);
            grid_cells[r_counter][ph_counter][th_counter]->setThID(th_counter);
            grid_cells[r_counter][ph_counter][th_counter]->resize(data_len + tmp_data_offset);
            tmp_cell = grid_cells[r_counter][ph_counter][th_counter];
            tmp_cell->setID(line_counter);
        }

        for(uint i = 0; i < data_offset; i++)
        {
            double tmp_data1 = 0;
            bin_reader.read((char *)&tmp_data1, 8);
            // cout << tmp_data1 << " ";
            tmp_cell->setData(i, tmp_data1);
        }

        updateVelocity(tmp_cell, param);

        if(uint(tmp_cell->getData(data_pos_id)) < 0 ||
           uint(tmp_cell->getData(data_pos_id)) > param.getMaxDustComponentChoice())
        {
            cout << "\nERROR: Cannot write to:\n Dust ID in grid exceeds maximum number "
                    "of dust choices available! "
                 << endl;
            return false;
        }

        updateDataRange(tmp_cell, param);

        double tmp_vol = getVolume(tmp_cell);
        total_gas_mass += getGasMassDensity(tmp_cell) * tmp_vol;
        cell_volume += tmp_vol;
        th_counter++;
    }

    bin_reader.close();

    if(max_cells != uint(line_counter))
    {
        cout << "\nERROR: Cannot write to:\n Number of read in cells do not match the "
                "maximal number of expected cells!"
             << endl;
        return false;
    }

    data_len += tmp_data_offset;
    data_offset += tmp_data_offset;

    max_len = 2 * Rmax;
    // min_len = listR[1] - listR[0];

    // cout << CLR_LINE;
    // cout << "- Loading spherical grid file          : done" << endl;

    return true;
}

bool CGridSpherical::writeGNUPlotFiles(string path, parameters & param)
{
    nrOfGnuPoints = param.getNrOfGnuPoints();
    nrOfGnuVectors = param.getNrOfGnuVectors();
    maxGridLines = param.getmaxGridLines();

    if(nrOfGnuPoints + nrOfGnuVectors == 0)
        return true;

    if(max_cells == 0)
    {
        cout << "\nERROR: Cannot plot spherical grid to Gnuplot file to:" << endl;
        cout << path;
        cout << "Not enough tree cells available! " << endl;
        return false;
    }

    plt_gas_dens = (!data_pos_gd_list.empty());  // 1
    plt_dust_dens = false;                       // param.getPlot(plIDnd) && (!data_pos_dd_list.empty()); // 2
    plt_gas_temp = (data_pos_tg != MAX_UINT);    // 3
    plt_dust_temp = (!data_pos_dt_list.empty()); // 4
    plt_rat = (!data_pos_aalg_list.empty());     // 5
    plt_delta = false;                           // param.getPlot(plIDdelta) && (data_pos_tg != MAX_UINT) &&
                       // (data_pos_mx != MAX_UINT) && (data_pos_td != MAX_UINT); // 6
    //plt_larm = false; // param.getPlot(plIDlarm) && (data_pos_tg != MAX_UINT) &&
    plt_mach = false; // param.getPlot(plIDmach) && (data_pos_vx != MAX_UINT) &&
                      // (data_pos_tg != MAX_UINT); // 7
    plt_disr = (!data_pos_adisr_list.empty());     // 8
    plt_max_disr = (!data_pos_max_adisr_list.empty());     // 9
    plt_param_modif = (!data_pos_param_modif_list.empty());     // 10
    plt_barnet_low_lower = (!data_pos_barnet_low_J_lower_list.empty());     // 11
    plt_barnet_low_upper = (!data_pos_barnet_low_J_upper_list.empty());     // 12
    plt_barnet_high_lower = (!data_pos_barnet_high_J_lower_list.empty());     // 13
    plt_barnet_high_upper = (!data_pos_barnet_high_J_upper_list.empty());     // 14
    plt_dg_lower = (!data_pos_dg_lower_list.empty());     // 15
    plt_dg_upper = (!data_pos_dg_upper_list.empty());     // 16
    plt_dg_10_lower = (!data_pos_dg_10_lower_list.empty());     // 17
    plt_dg_10_upper = (!data_pos_dg_10_upper_list.empty());     // 18
    plt_abs_ini = (!data_pos_abs_ini_list.empty()); // 19
    plt_amaxJB_Lar = (!data_pos_amaxJB_Lar_list.empty());     // 20
 

    plt_mag = (data_pos_mx != MAX_UINT); // 0
    plt_vel = (data_pos_vx != MAX_UINT); // 1

    if(nrOfGnuPoints <= 1)
    {
        nrOfGnuPoints = max_cells / 10;

        plt_gas_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_rat = false;
        plt_delta = false;
        //plt_larm = false;
        plt_mach = false;
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
        plt_abs_ini = false;
        plt_amaxJB_Lar = false;
    }
    else
        nrOfGnuPoints = max_cells / nrOfGnuPoints;

    if(nrOfGnuVectors <= 1)
    {
        nrOfGnuVectors = max_cells / 10;
        plt_mag = false;
        plt_vel = false;
    }
    else
        nrOfGnuVectors = max_cells / nrOfGnuVectors;

    if(nrOfGnuPoints == 0)
        nrOfGnuPoints = 1;

    if(nrOfGnuVectors == 0)
        nrOfGnuVectors = 1;

    stringstream point_header, vec_header, basic_grid_l0, basic_grid_l1;

    string grid_filename = path + "grid_geometry.plt";
    string dens_gas_filename = path + "grid_gas_density.plt";
    string dens_dust_filename = path + "grid_dust_density.plt";
    string temp_gas_filename = path + "grid_gas_temp.plt";
    string temp_dust_filename = path + "grid_dust_temp.plt";
    string abs_ini_filename = path + "grid_abs_ini.py";
    string rat_filename = path + "grid_RAT.plt";
    string delta_filename = path + "grid_data.dat";
    //string larm_filename = path + "grid_mag.plt";
    string mach_filename = path + "grid_vel.plt";
    string ratd_filename = path + "grid_RATD.py";
    string max_ratd_filename = path + "grid_max_RATD.py";
    string param_modif_filename = path + "grid_param_modif.py";
    string barnet_low_lower_filename = path + "grid_amin_aJ_lowJ.py";
    string barnet_low_upper_filename = path + "grid_amax_aJ_lowJ.py";
    string barnet_high_lower_filename = path + "grid_amin_aJ_highJ.py";
    string barnet_high_upper_filename = path + "grid_amax_aJ_highJ.py";
    string dg_lower_filename = path + "grid_amin_dg_0.5.py";
    string dg_upper_filename = path + "grid_amax_dg_1.py";
    string dg_10_lower_filename = path + "grid_amin_dg_1.py";
    string dg_10_upper_filename = path + "grid_amax_dg_1.py";
    string amaxJB_Lar_filename = path + "grid_amaxJB_Lar.plt";
    string mag_filename = path + "grid_mag.py";
    string vel_filename = path + "grid_vel.py";


    ofstream point_fields[21];

    ofstream vec_fields[2];

    point_fields[0].open(grid_filename.c_str());

    if(point_fields[0].fail())
    {
        cout << "\nERROR: Cannot write to:\n" << grid_filename << endl;
        return false;
    }

    if(plt_gas_dens)
    {
        point_fields[1].open(dens_gas_filename.c_str());

        if(point_fields[1].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << dens_gas_filename << endl;
            return false;
        }
    }

    if(plt_dust_dens)
    {
        point_fields[2].open(dens_dust_filename.c_str());

        if(point_fields[2].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << dens_dust_filename << endl;
            return false;
        }
    }

    if(plt_gas_temp)
    {
        point_fields[3].open(temp_gas_filename.c_str());

        if(point_fields[3].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << temp_gas_filename << endl;
            return false;
        }
    }

    if(plt_dust_temp)
    {
        point_fields[4].open(temp_dust_filename.c_str());

        if(point_fields[4].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << temp_dust_filename << endl;
            return false;
        }
    }

    if(plt_rat)
    {
        point_fields[5].open(rat_filename.c_str());

        if(point_fields[5].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << rat_filename << endl;
            return false;
        }
    }

    if(plt_delta)
    {
        point_fields[6].open(delta_filename.c_str());

        if(point_fields[6].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << delta_filename << endl;
            return false;
        }
    }

    //if(plt_larm)
    //{
    //    point_fields[7].open(larm_filename.c_str());

    //    if(point_fields[7].fail())
    //    {
    //        cout << "\nERROR: Cannot write to:\n" << larm_filename << endl;
    //        return false;
    //    }
    //}
    

    if(plt_mach)
    {
        point_fields[7].open(mach_filename.c_str());

        if(point_fields[7].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << mach_filename << endl;
            return false;
        }
    }

    if(plt_disr)
    {
        point_fields[8].open(ratd_filename.c_str());

        if(point_fields[8].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << ratd_filename << endl;
            return false;
        }
    }

    if(plt_max_disr)
    {
        point_fields[9].open(max_ratd_filename.c_str());

        if(point_fields[9].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << max_ratd_filename << endl;
            return false;
        }
    }

    if(plt_param_modif)
    {
        point_fields[10].open(param_modif_filename.c_str());

        if(point_fields[10].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << param_modif_filename << endl;
            return false;
        }
    }
    
    if(plt_barnet_low_lower)
    {
        point_fields[11].open(barnet_low_lower_filename.c_str());

        if(point_fields[11].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << barnet_low_lower_filename << endl;
            return false;
        }
    }
    
    if(plt_barnet_low_upper)
    {
        point_fields[12].open(barnet_low_upper_filename.c_str());

        if(point_fields[12].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << barnet_low_upper_filename << endl;
            return false;
        }
    }

    
    if(plt_barnet_high_lower)
    {
        point_fields[13].open(barnet_high_lower_filename.c_str());

        if(point_fields[13].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << barnet_high_lower_filename << endl;
            return false;
        }
    }
    
    if(plt_barnet_high_upper)
    {
        point_fields[14].open(barnet_high_upper_filename.c_str());

        if(point_fields[14].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << barnet_high_upper_filename << endl;
            return false;
        }
    }
    
    if(plt_dg_lower)
    {
        point_fields[15].open(dg_lower_filename.c_str());

        if(point_fields[15].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << dg_lower_filename << endl;
            return false;
        }
    }
    
    if(plt_dg_upper)
    {
        point_fields[16].open(dg_upper_filename.c_str());

        if(point_fields[16].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << dg_upper_filename << endl;
            return false;
        }
    }
    
    if(plt_dg_10_lower)
    {
        point_fields[17].open(dg_10_lower_filename.c_str());

        if(point_fields[17].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << dg_10_lower_filename << endl;
            return false;
        }
    }
    
    if(plt_dg_10_upper)
    {
        point_fields[18].open(dg_10_upper_filename.c_str());

        if(point_fields[18].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << dg_10_upper_filename << endl;
            return false;
        }
    }
    
    
    if(plt_abs_ini)
    {
        point_fields[19].open(abs_ini_filename.c_str());

        if(point_fields[19].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << abs_ini_filename << endl;
            return false;
        }
    }
    
    if(plt_amaxJB_Lar)
    {
        point_fields[20].open(amaxJB_Lar_filename.c_str());

        if(point_fields[20].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << amaxJB_Lar_filename << endl;
            return false;
        }
    }
    
    if(plt_mag)
    {
        vec_fields[0].open(mag_filename.c_str());

        if(vec_fields[0].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << mag_filename << endl;
            return false;
        }
    }

    if(plt_vel)
    {
        vec_fields[1].open(vel_filename.c_str());

        if(vec_fields[1].fail())
        {
            cout << "\nERROR: Cannot write to:\n" << vel_filename << endl;
            return false;
        }
    }

    line_counter = 0;
    char_counter = 0;

    // Grid boundaries
    for(uint i_th = 0; i_th < N_th; i_th++)
    {
        double th = listTh[i_th];
        double Nstep = 45.0;

        for(uint i_ph = 0; i_ph < Nstep; i_ph++)
        {
            double tmp_dph = PIx2 / Nstep;
            Vector3D p1, p2, dist;

            // inner sphere
            p1 = Vector3D(Rmin, i_ph * tmp_dph, th);
            p2 = Vector3D(Rmin, (i_ph + 1) * tmp_dph, th);

            p1.spher2cart();
            p2.spher2cart();
            dist = p2 - p1;

            basic_grid_l0 << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " "
                          << float(dist.X()) << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;

            // outer sphere
            p1 = Vector3D(Rmax, i_ph * tmp_dph, th);
            p2 = Vector3D(Rmax, (i_ph + 1) * tmp_dph, th);

            p1.spher2cart();
            p2.spher2cart();
            dist = p2 - p1;

            basic_grid_l1 << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " "
                          << float(dist.X()) << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
        }
    }

    for(uint i_ph = 0; i_ph < N_ph; i_ph++)
    {
        double ph = listPh[i_ph];
        double Nstep = 25.0;

        for(uint i_th = 0; i_th < Nstep; i_th++)
        {
            double tmp_dth = PI / Nstep;
            Vector3D p1, p2, dist;

            // inner sphere
            p1 = Vector3D(Rmin, ph, i_th * tmp_dth);
            p2 = Vector3D(Rmin, ph, (i_th + 1) * tmp_dth);

            p1.spher2cart();
            p2.spher2cart();
            dist = p2 - p1;

            basic_grid_l0 << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " "
                          << float(dist.X()) << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;

            // outer sphere
            p1 = Vector3D(Rmax, ph, i_th * tmp_dth);
            p2 = Vector3D(Rmax, ph, (i_th + 1) * tmp_dth);

            p1.spher2cart();
            p2.spher2cart();
            dist = p2 - p1;

            basic_grid_l1 << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " "
                          << float(dist.X()) << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
        }
    }

    cout << CLR_LINE;

    point_header.str("");
    point_header << "reset" << endl;
    point_header << "#set terminal postscript" << endl;
    point_header << "#set output \'filename.plt\'" << endl;
    point_header << "set ticslevel 0" << endl;
    point_header << "set size ratio -1" << endl;
    point_header << "set view 45,45" << endl;

    point_header << "set xlabel \'x[m]\'" << endl;
    point_header << "set ylabel \'y[m]\'" << endl;
    point_header << "set zlabel \'z[m]\'" << endl;

    point_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    point_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    point_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    point_header << "set style arrow 1 nohead lt 2 lc rgb 0x0000ff lw 0.5" << endl;
    point_header << "set style arrow 2 nohead lt 2 lc rgb 0x0054dd lw 0.5" << endl;

    point_header << "set style line 1 pt 1 ps variable lt palette" << endl;

    point_header << "set grid" << endl;
    point_header << "set nokey" << endl;

    // 0 spherical grid
    point_fields[0] << point_header.str();
    point_fields[0] << "set title \'3D spherical grid geometry\' font \'Arial,12\'" << endl;
    point_fields[0] << "set style arrow 3 nohead ls 1 lw 0.5 lc rgb 0x550066" << endl;
    point_fields[0] << "splot '-' with vectors as 3,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 1 gas density
    point_fields[1] << point_header.str();
    point_fields[1] << "set title \'3D gas number density distribution (min: " << min_gas_dens
                    << "[m^-3]; max: " << max_gas_dens << "[m^-3])\' font \'Arial,12\'" << endl;
    point_fields[1] << "set cblabel \'gas density[m^-3]\'" << endl;
    point_fields[1] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1, 3 1 1 0)" << endl;

    if(min_gas_dens == 0 && max_gas_dens == 0)
    {
        min_gas_dens = 0.1;
        max_gas_dens = 1;
    }

    if(min_gas_dens == 0)
    {
        min_gas_dens = 0.001 * max_gas_dens;
    }

    if(min_gas_dens / max_gas_dens > 0.9)
        min_gas_dens = 0.9 * max_gas_dens;

    point_fields[1] << "set cbrange[" << log10(min_gas_dens) << ":" << log10(max_gas_dens) << "]" << endl;
    point_fields[1] << "set format cb \'%.02g\'" << endl;

    point_fields[1] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 2 dust density
    point_fields[2] << point_header.str();
    point_fields[2] << "set title \'3D gas number density distribution (min: " << min_dust_dens
                    << "[m^-3]; max: " << max_dust_dens << "[m^-3])\' font \'Arial,12\'" << endl;
    point_fields[2] << "set cblabel \'gas density[m^-3]\'" << endl;
    point_fields[2] << "set palette defined (0 0.5 0 0, 1 0 0 1, 2 0 1 1)" << endl;

    if(min_dust_dens == 0 && max_dust_dens == 0)
    {
        min_dust_dens = 0.1;
        max_dust_dens = 1;
    }

    if(min_dust_dens == 0)
    {
        min_dust_dens = 0.001 * max_dust_dens;
    }

    if(min_dust_dens / max_dust_dens > 0.9)
        min_dust_dens = 0.9 * max_dust_dens;

    point_fields[2] << "set cbrange[" << log10(min_dust_dens) << ":" << log10(max_dust_dens) << "]" << endl;
    point_fields[2] << "set format cb \'%.02g\'" << endl;

    point_fields[2] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 3 gas_temp
    point_fields[3] << point_header.str();
    point_fields[3] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[3] << "set title \'3D gas temperature distribution (min: " << min_gas_temp
                    << "[K]; max: " << max_gas_temp << "[K])\' font \'Arial,12\'" << endl;
    point_fields[3] << "set cblabel \'temperature [K]\'" << endl;

    if(min_gas_temp == 0 && max_gas_temp == 0)
    {
        min_gas_temp = 0.1;
        max_gas_temp = 1;
    }

    if(min_gas_temp / max_gas_temp > 0.90)
        min_gas_temp = 0.9 * max_gas_temp;

    point_fields[3] << "set cbrange[" << float(min_gas_temp) << ":" << float(max_gas_temp) << "]" << endl;
    point_fields[3] << "set format cb \'%.03g\'" << endl;

    point_fields[3] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 4 dust temp
    point_fields[4] << point_header.str();
    point_fields[4] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[4] << "set title \'3D dust temperature distribution (min: " << min_dust_temp
                    << "[K]; max: " << max_dust_temp << "[K])\' font \'Arial,12\'" << endl;
    point_fields[4] << "set cblabel \'temperature [K]\'" << endl;

    if(min_dust_temp == 0 && max_dust_temp == 0)
    {
        min_dust_temp = 0.1;
        max_dust_temp = 1;
    }

    if(min_dust_temp / max_dust_temp > 0.9)
        min_dust_temp = 0.9 * max_dust_temp;

    point_fields[4] << "set cbrange[" << float(min_dust_temp) << ":" << float(max_dust_temp) << "]" << endl;
    point_fields[4] << "set format cb \'%.03g\'" << endl;

    point_fields[4] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    // 5 rat
    point_fields[5] << point_header.str();
    point_fields[5] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[5] << "set title \'3D aligned radii distribution (min ID: " << aalg_min
                    << "; max ID: " << aalg_max << ")\' font \'Arial,12\'" << endl;

    point_fields[5] << "set cblabel \'aligned radius ID\'" << endl;

    if(aalg_min == aalg_max)
        aalg_max = 1.01 * aalg_min;

    point_fields[5] << "set cbrange[" << aalg_min << ":" << aalg_max << "]" << endl;

    point_fields[5] << "set format cb \'%.03g\'" << endl;

    point_fields[5] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

// 9 ratd
    point_fields[9] << point_header.str();
    point_fields[9] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[9] << "set title \'3D disrupted radii distribution (min ID: " << adisr_min
                    << "; max ID: " << adisr_max << ")\' font \'Arial,12\'" << endl;

    point_fields[9] << "set cblabel \'disrupted radius ID\'" << endl;

    if(adisr_min == adisr_max)
        adisr_max = 1.01 * adisr_min;

    point_fields[9] << "set cbrange[" << adisr_min << ":" << adisr_max << "]" << endl;

    point_fields[9] << "set format cb \'%.03g\'" << endl;

    point_fields[9] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;


// 10 ratd_max
    point_fields[10] << point_header.str();
    point_fields[10] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[10] << "set title \'3D maximum disrupted radii distribution (min ID: " << max_adisr_min
                    << "; max ID: " << max_adisr_max << ")\' font \'Arial,12\'" << endl;

    point_fields[10] << "set cblabel \'disrupted radius ID\'" << endl;

    if(max_adisr_min == max_adisr_max)
        max_adisr_max = 1.01 * max_adisr_min;

    point_fields[10] << "set cbrange[" << max_adisr_min << ":" << max_adisr_max << "]" << endl;

    point_fields[10] << "set format cb \'%.03g\'" << endl;

    point_fields[10] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

// 11 size_param_modify
    point_fields[11] << point_header.str();
    point_fields[11] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[11] << "set title \'3D size distribution slope (min ID: " << size_param_modif_min
                    << "; max ID: " << size_param_modif_max << ")\' font \'Arial,12\'" << endl;

    point_fields[11] << "set cblabel \'power-law size distribution ID\'" << endl;

    if(size_param_modif_min == size_param_modif_max)
        size_param_modif_max = 1.01 * size_param_modif_min;

    point_fields[11] << "set cbrange[" << size_param_modif_min << ":" << size_param_modif_max << "]" << endl;

    point_fields[11] << "set format cb \'%.03g\'" << endl;

    point_fields[11] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;
 
 
// 12 abarnet_low_J_lower
    point_fields[12] << point_header.str();
    point_fields[12] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[12] << "set title \'3D barnet radii at low J distribution [lower limit] (min ID: " << abar_low_lower_min
                    << "; max ID: " << abar_low_lower_max << ")\' font \'Arial,12\'" << endl;

    point_fields[12] << "set cblabel \'lower barnett radius at low J ID\'" << endl;

    if(abar_low_lower_min == abar_low_lower_max)
        abar_low_lower_max = 1.01 * abar_low_lower_min;

    point_fields[12] << "set cbrange[" << abar_low_lower_min << ":" << abar_low_lower_max << "]" << endl;

    point_fields[12] << "set format cb \'%.03g\'" << endl;

    point_fields[12] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

// 13 abarnet_low_J_upper
    point_fields[13] << point_header.str();
    point_fields[13] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[13] << "set title \'3D barnet radii at low J distribution [upper limit] (min ID: " << abar_low_upper_min
                    << "; max ID: " << abar_low_upper_max << ")\' font \'Arial,12\'" << endl;

    point_fields[13] << "set cblabel \'upper barnett radius at low J ID\'" << endl;

    if(abar_low_upper_min == abar_low_upper_max)
        abar_low_upper_max = 1.01 * abar_low_upper_min;

    point_fields[13] << "set cbrange[" << abar_low_upper_min << ":" << abar_low_upper_max << "]" << endl;

    point_fields[13] << "set format cb \'%.03g\'" << endl;

    point_fields[13] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;
 
 
// 14 abarnet_high_J_lower
    point_fields[14] << point_header.str();
    point_fields[14] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[14] << "set title \'3D barnet radii at high J distribution [lower limit] (min ID: " << abar_high_lower_min
                    << "; max ID: " << abar_high_lower_max << ")\' font \'Arial,12\'" << endl;

    point_fields[14] << "set cblabel \'lower barnett radius at low J ID\'" << endl;

    if(abar_high_lower_min == abar_high_lower_max)
        abar_high_lower_max = 1.01 * abar_high_lower_min;

    point_fields[14] << "set cbrange[" << abar_high_lower_min << ":" << abar_high_lower_max << "]" << endl;

    point_fields[14] << "set format cb \'%.03g\'" << endl;

    point_fields[14] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

// 15 abarnet_high_J_upper
    point_fields[15] << point_header.str();
    point_fields[15] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[15] << "set title \'3D barnet radii at high J distribution [upper limit] (min ID: " << abar_high_upper_min
                    << "; max ID: " << abar_high_upper_max << ")\' font \'Arial,12\'" << endl;

    point_fields[15] << "set cblabel \'upper barnett radius at low J ID\'" << endl;

    if(abar_high_upper_min == abar_high_upper_max)
        abar_high_upper_max = 1.01 * abar_high_upper_min;

    point_fields[15] << "set cbrange[" << abar_high_upper_min << ":" << abar_high_upper_max << "]" << endl;

    point_fields[15] << "set format cb \'%.03g\'" << endl;

    point_fields[15] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;
    
// 16 adg_lower
    point_fields[16] << point_header.str();
    point_fields[16] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[16] << "set title \'3D threshold of the lower limit of tau_mag < tau_gas (min ID: " << adg_lower_min
                    << "; max ID: " << adg_lower_max << ")\' font \'Arial,12\'" << endl;

    point_fields[16] << "set cblabel \'lower threshold of the lower limit of tau_mag < tau_gas ID\'" << endl;

    if(adg_lower_min == adg_lower_max)
        adg_lower_max = 1.01 * adg_lower_min;

    point_fields[16] << "set cbrange[" << adg_lower_min << ":" << adg_lower_max << "]" << endl;

    point_fields[16] << "set format cb \'%.03g\'" << endl;

    point_fields[16] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

// 17 adg_upper
    point_fields[17] << point_header.str();
    point_fields[17] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[17] << "set title \'3D threshold of the upper limit of tau_mag < tau_gas (min ID: " << adg_upper_min
                    << "; max ID: " << adg_upper_max << ")\' font \'Arial,12\'" << endl;

    point_fields[17] << "set cblabel \'upper threshold of the upper limit of tau_mag < tau_gas ID\'" << endl;

    if(adg_upper_min == adg_upper_max)
        adg_upper_max = 1.01 * adg_upper_min;

    point_fields[17] << "set cbrange[" << adg_upper_min << ":" << adg_upper_max << "]" << endl;

    point_fields[17] << "set format cb \'%.03g\'" << endl;

    point_fields[17] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;
    
// 18 adg_10_lower
    point_fields[18] << point_header.str();
    point_fields[18] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[18] << "set title \'3D threshold of the lower limit of tau_mag < tau_gas/10 (min ID: " << adg_10_lower_min
                    << "; max ID: " << adg_10_lower_max << ")\' font \'Arial,12\'" << endl;

    point_fields[18] << "set cblabel \'lower threshold of the lower limit of tau_mag < tau_gas/10 ID\'" << endl;

    if(adg_10_lower_min == adg_10_lower_max)
        adg_10_lower_max = 1.01 * adg_10_lower_min;

    point_fields[18] << "set cbrange[" << adg_10_lower_min << ":" << adg_10_lower_max << "]" << endl;

    point_fields[18] << "set format cb \'%.03g\'" << endl;

    point_fields[18] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;

// 19 adg_10_upper
    point_fields[19] << point_header.str();
    point_fields[19] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[19] << "set title \'3D threshold of the upper limit of tau_mag < tau_gas/10 (min ID: " << adg_10_upper_min
                    << "; max ID: " << adg_10_upper_max << ")\' font \'Arial,12\'" << endl;

    point_fields[19] << "set cblabel \'upper threshold for f_highJ = 0.5 at low J ID\'" << endl;

    if(adg_10_upper_min == adg_10_upper_max)
        adg_10_upper_max = 1.01 * adg_10_upper_min;

    point_fields[19] << "set cbrange[" << adg_10_upper_min << ":" << adg_10_upper_max << "]" << endl;

    point_fields[19] << "set format cb \'%.03g\'" << endl;

    point_fields[19] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;

    vec_header.str("");
    vec_header << "reset" << endl;
    vec_header << "#set terminal postscript" << endl;
    vec_header << "#set output \'\'" << endl;
    vec_header << "set ticslevel 0" << endl;
    vec_header << "set size ratio -1" << endl;
    vec_header << "set view 45,45" << endl;

    vec_header << "set xlabel \'x[m]\'" << endl;
    vec_header << "set ylabel \'y[m]\'" << endl;
    vec_header << "set zlabel \'z[m]\'" << endl;

    vec_header << "set xrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set yrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;
    vec_header << "set zrange[" << -1.01 * Rmax << ":" << 1.01 * Rmax << "]" << endl;

    vec_header << "set style arrow 1 nohead ls 1 lw 1 lc rgb 0x0000cc" << endl;
    vec_header << "set style arrow 2 nohead ls 1 lw 1 lc rgb 0x5500dd" << endl;
    vec_header << "set style arrow 3 ls 1 lw 1 lc palette" << endl;

    vec_header << "set grid" << endl;
    vec_header << "set nokey" << endl;
    
    
// 20 abs_ini: initial absorption rate of dust grains inside the grid
    point_fields[20] << point_header.str();
    point_fields[20] << "set palette defined (0 0.05 0 0, 0.4 1 0 0, 0.7 1 1 0, 1 1 1 0.5)" << endl;

    point_fields[20] << "set title \'3D initial dust absorption rate (min: " << min_abs_ini
                    << "; max: " << max_abs_ini << ")\' font \'Arial,12\'" << endl;
    point_fields[20] << "set cblabel \'absorption rate\'" << endl;

    if(min_abs_ini == 0 && max_abs_ini == 0)
    {
        min_abs_ini = 0.1;
        max_abs_ini = 1;
    }

    if(min_abs_ini / max_abs_ini > 0.9)
        min_abs_ini = 0.9 * max_abs_ini;

    point_fields[20] << "set cbrange[" << float(min_abs_ini) << ":" << float(max_abs_ini) << "]" << endl;
    point_fields[20] << "set format cb \'%.03g\'" << endl;

    point_fields[20] << "splot  '-' w p ls 1,'-' with vectors as 2,'-' with vectors as 1" << endl;
    
    
// 0 mag
    vec_fields[0] << vec_header.str();
    vec_fields[0] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 0 0.1 0.9)" << endl;

    if(min_mag == 0 && max_mag == 0)
    {
        min_mag = 1e-45;
        max_mag = 2e-45;
    }

    vec_fields[0] << "set title \'3D mag. field distribution (min:" << log10(min_mag)
                  << " log10([T]); max:" << log10(max_mag) << " log10([T])  \' font \'Arial,12\'" << endl;

    if(min_mag / max_mag > 0.9)
        min_mag = 0.9 * max_mag;

    vec_fields[0] << "set cbrange[" << log10(min_mag) << ":" << log10(max_mag) << "]" << endl;
    vec_fields[0] << "set format cb \'%.02g\'" << endl;
    vec_fields[0] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' "
                     "with vectors as 1"
                  << endl;

// 1 vel
    vec_fields[1] << vec_header.str();
    vec_fields[1] << "set palette defined (0 1 0 0, 0.5 0.0 0.9 0,  0.75 0.0 0.9 1, 0.9 0 0.1 0.9)" << endl;

    if(min_vel == 0 && max_vel == 0)
    {
        min_vel = 1e-45;
        max_vel = 1e-45;
    }

    vec_fields[1] << "set title \'3D vel. field directions (min:" << log10(min_vel)
                  << " log10(m/s); max:" << log10(max_vel) << " log10(m/s)\' font \'Arial,12\'" << endl;

    if(min_vel / max_vel > 0.9)
        min_vel = 0.9 * max_vel;

    vec_fields[1] << "set cbrange[" << float(log10(min_vel)) << ":" << float(log10(max_vel)) << "]" << endl;
    vec_fields[1] << "set format cb \'%.03g\'" << endl;
    vec_fields[1] << "splot  \'-\' with vectors as 3, \'-\' with vectors as 2, \'-\' "
                     "with vectors as 1"
                  << endl;

    if(plt_gas_dens)
        point_fields[0] << "0 0 0 " << 0.5 << " " << getGasDensity(center_cell) << endl;

    line_counter = 0;

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        double size = 0.5 + 1.5 * double(i_r) / double(N_r);

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_th = 0; i_th < N_th; i_th++)
            {
                cell_sp * tmp_cell_pos = grid_cells[i_r][i_ph][i_th];

                Vector3D c = getCenter(tmp_cell_pos);

                line_counter++;

                double scale = 0;

                double dr = listR[i_r + 1] - listR[i_r];
                double dph = listPh[i_ph + 1] - listPh[i_ph];
                double dth = listTh[i_th + 1] - listTh[i_th];

                if(line_counter % nrOfGnuPoints == 0)
                {
                    if(plt_gas_dens)
                    {
                        double dens = getGasDensity(tmp_cell_pos);

                        if(dens > 0)
                            point_fields[1] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size)
                                            << " " << log10(dens) << endl;
                    }

                    if(plt_gas_temp)
                    {
                        double Tg = getGasTemperature(tmp_cell_pos);
                        point_fields[3] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << Tg << endl;
                    }

                    if(plt_gas_temp)
                    {
                        double Td = getDustTemperature(tmp_cell_pos);
                        point_fields[4] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << Td << endl;
                    }

                    if(plt_rat)
                    {
                        double a_alg = getAlignedRadius(tmp_cell_pos, 0);
                        point_fields[5] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_alg << endl;
                    }

                    if(plt_disr)
                    {
                        double a_disr = getDisruptRadius(tmp_cell_pos, 0);
                        point_fields[9] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_disr << endl;
                    }

                    if(plt_max_disr)
                    {
                        double max_a_disr = getMaxDisruptRadius(tmp_cell_pos, 0);
                        point_fields[10] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << max_a_disr << endl;
                    }

                    if(plt_param_modif)
                    {
                        double size_param_modify = getSizeParamModify(tmp_cell_pos, 0);
                        point_fields[11] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << size_param_modify << endl;
                    }
                    
                    if(plt_barnet_low_lower)
                    {
                        double a_bar_low_lower = getBarnetLowLowerRadius(tmp_cell_pos, 0);
                        point_fields[12] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_bar_low_lower << endl;
                    }
               
               		if(plt_barnet_low_upper)
                    {
                        double a_bar_low_upper = getBarnetLowUpperRadius(tmp_cell_pos, 0);
                        point_fields[13] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_bar_low_upper << endl;
                    }
                    
                    if(plt_barnet_high_lower)
                    {
                        double a_bar_high_lower = getBarnetHighLowerRadius(tmp_cell_pos, 0);
                        point_fields[14] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_bar_high_lower << endl;
                    }
               
               		if(plt_barnet_high_upper)
                    {
                        double a_bar_high_upper = getBarnetHighUpperRadius(tmp_cell_pos, 0);
                        point_fields[15] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_bar_high_upper << endl;
                    }
                    
                    if(plt_dg_lower)
                    {
                        double a_dg_lower = getDGLowerRadius(tmp_cell_pos, 0);
                        point_fields[16] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_dg_lower << endl;
                    }
               
               		if(plt_dg_upper)
                    {
                        double a_dg_upper = getDGUpperRadius(tmp_cell_pos, 0);
                        point_fields[17] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_dg_upper << endl;
                    }
                    
                    
                    if(plt_dg_10_lower)
                    {
                        double a_dg_10_lower = getDG10LowerRadius(tmp_cell_pos, 0);
                        point_fields[18] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_dg_10_lower << endl;
                    }
               
               		if(plt_dg_10_upper)
                    {
                        double a_dg_10_upper = getDG10UpperRadius(tmp_cell_pos, 0);
                        point_fields[19] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << a_dg_10_upper << endl;
                    }
                    
                    
                    if(plt_abs_ini)
                    {
                        double abs_ini = getQBOffset(tmp_cell_pos);
                        point_fields[20] << c.X() << " " << c.Y() << " " << c.Z() << " " << float(size) << " "
                                        << abs_ini << endl;
                    }
                }

                if(line_counter % nrOfGnuVectors == 0)
                {
                    if(plt_mag)
                    {
                        double mx = getMagField(tmp_cell_pos).X();
                        double my = getMagField(tmp_cell_pos).Y();
                        double mz = getMagField(tmp_cell_pos).Z();

                        double b_len = sqrt(mx * mx + my * my + mz * mz);

                        if(b_len != 0)
                        {
                            mx = 0.45 * min_len * size * mx / b_len;
                            my = 0.45 * min_len * size * my / b_len;
                            mz = 0.45 * min_len * size * mz / b_len;

                            vec_fields[0] << float(c.X() - mx) << " " << float(c.Y() - my) << " "
                                          << float(c.Z() - mz) << " " << float(2.0 * mx) << " "
                                          << float(2.0 * my) << " " << float(2.0 * mz) << " "
                                          << float(log10(b_len)) << endl;
                        }
                    }

                    if(plt_vel)
                    {
                        double vx = getVelocityField(tmp_cell_pos).X();
                        double vy = getVelocityField(tmp_cell_pos).Y();
                        double vz = getVelocityField(tmp_cell_pos).Z();

                        double v_len = sqrt(vx * vx + vy * vy + vz * vz);

                        if(v_len != 0)
                        {
                            vx = 0.45 * min_len * size * vx / v_len;
                            vy = 0.45 * min_len * size * vy / v_len;
                            vz = 0.45 * min_len * size * vz / v_len;

                            vec_fields[1] << float(c.X() - vx) << " " << float(c.Y() - vy) << " "
                                          << float(c.Z() - vz) << " " << float(2.0 * vx) << " "
                                          << float(2.0 * vy) << " " << float(2.0 * vz) << " "
                                          << float(log10(v_len)) << endl;
                        }
                    }
                }
            }
        }
    }

    for(uint pos = 1; pos < 6; pos++)
    {
        point_fields[pos] << "\ne\n"
                          << basic_grid_l0.str() << "\ne\n"
                          << basic_grid_l1.str() << "\ne" << endl;
    }

    point_fields[9] << "\ne\n" //GRIDadisr
                      << basic_grid_l0.str() << "\ne\n"
                      << basic_grid_l1.str() << "\ne" << endl;

    point_fields[10] << "\ne\n" //GRIDadisr_max
                      << basic_grid_l0.str() << "\ne\n"
                      << basic_grid_l1.str() << "\ne" << endl;

    point_fields[11] << "\ne\n" //GRIDparam_modif
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                    
    point_fields[12] << "\ne\n" //GRIDabar_low_lower
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
	point_fields[13] << "\ne\n" //GRIDabar_low_upper
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;

    point_fields[14] << "\ne\n" //GRIDabar_high_lower
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
	point_fields[15] << "\ne\n" //GRIDabar_high_upper
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
    point_fields[16] << "\ne\n" //GRIDadg_lower
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
	point_fields[17] << "\ne\n" //GRIDadg_upper
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
    point_fields[18] << "\ne\n" //GRIDadg_10_lower
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
	point_fields[19] << "\ne\n" //GRIDadg_10_upper
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
	point_fields[20] << "\ne\n" //GRIDabs_ini
                      << basic_grid_l0.str() << "\ne\n"
                     << basic_grid_l1.str() << "\ne" << endl;
                     
    for(uint pos = 0; pos < 2; pos++)
    {
        vec_fields[pos] << "\ne\n" << basic_grid_l0.str() << "\ne\n" << basic_grid_l1.str() << "\ne" << endl;
    }

    stringstream buffer;
    buffer.str("");

    for(uint i_th = 0; i_th <= N_th; i_th++)
    {
        double th = listTh[i_th];

        for(uint i_r = 0; i_r <= N_r; i_r++)
        {
            double r = listR[i_r];
            double Nstep = double(int(10 + 30.0 * double(i_r) / double(N_r) + 0.5));

            for(uint i_ph = 0; i_ph < Nstep; i_ph++)
            {
                double tmp_dph = PIx2 / Nstep;
                Vector3D p1(r, i_ph * tmp_dph, th);
                Vector3D p2(r, (i_ph + 1) * tmp_dph, th);

                p1.spher2cart();
                p2.spher2cart();
                Vector3D dist = p2 - p1;

                buffer << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " "
                       << float(dist.X()) << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
            }
        }
    }

    point_fields[0] << buffer.str() << "\ne" << endl;

    buffer.str("");

    for(uint i_ph = 0; i_ph <= N_ph; i_ph++)
    {
        double ph = listPh[i_ph];

        for(uint i_r = 0; i_r <= N_r; i_r++)
        {
            double r = listR[i_r];
            double Nstep = double(int(10 + 30.0 * double(i_r) / double(N_r) + 0.5));

            for(uint i_th = 0; i_th < Nstep; i_th++)
            {
                double tmp_dth = PI / Nstep;
                Vector3D p1(r, ph, i_th * tmp_dth);
                Vector3D p2(r, ph, (i_th + 1) * tmp_dth);

                p1.spher2cart();
                p2.spher2cart();
                Vector3D dist = p2 - p1;

                buffer << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " "
                       << float(dist.X()) << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
            }
        }
    }

    point_fields[0] << buffer.str() << "\ne" << endl;

    buffer.str("");

    for(uint i_ph = 0; i_ph <= N_ph; i_ph++)
    {
        double ph = listPh[i_ph];

        for(uint i_th = 0; i_th <= N_th; i_th++)
        {
            double th = listTh[i_th];
            Vector3D p1(Rmin, ph, th);
            Vector3D p2(Rmax, ph, th);

            p1.spher2cart();
            p2.spher2cart();
            Vector3D dist = p2 - p1;

            buffer << float(p1.X()) << " " << float(p1.Y()) << " " << float(p1.Z()) << " " << float(dist.X())
                   << " " << float(dist.Y()) << " " << float(dist.Z()) << endl;
        }
    }

    point_fields[0] << buffer.str() << "\ne" << endl;

    for(uint pos = 0; pos < 8; pos++)
        point_fields[pos].close();

    point_fields[9].close();
    point_fields[10].close();
    point_fields[11].close();
    point_fields[12].close();
    point_fields[13].close();
    point_fields[14].close();
    point_fields[15].close();
    point_fields[16].close();

    for(uint pos = 0; pos < 2; pos++)
        vec_fields[pos].close();

    cout << "- Writing of Gnuplot files             : done" << endl;
    return true;
}

bool CGridSpherical::saveBinaryGridFile(string filename, ushort id, ushort data_size)
{
    if(max_cells == 0)
    {
        cout << "\nERROR: Cannot save spherical grid file to:" << endl;
        cout << filename;
        cout << "Not enough cells available! " << endl;
        return false;
    }

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n Cannot open spherical grid file:" << endl;
        cout << filename;
        return false;
    }

    bin_writer.write((char *)&id, 2);
    bin_writer.write((char *)&data_size, 2);

    if(dataID == GRID_ID_SPH)
    {
        for(uint i = 0; i < data_offset; i++)
        {
            ushort tmp_ids = data_ids[i];
            bin_writer.write((char *)&tmp_ids, 2);
        }
    }
    else
    {
        cout << "\nERROR: Cannot save spherical grid file to:" << endl;
        cout << filename;
        cout << "A spherical grid requires an ID of " << GRID_ID_SPH << "!" << endl;
        return false;
    }

    bin_writer.write((char *)&Rmin, 8);
    bin_writer.write((char *)&Rmax, 8);
    bin_writer.write((char *)&N_r, 2);
    bin_writer.write((char *)&N_ph, 2);
    bin_writer.write((char *)&N_th, 2);
    bin_writer.write((char *)&log_factorR, 8);
    bin_writer.write((char *)&log_factorPh, 8);
    bin_writer.write((char *)&log_factorTh, 8);
    if(log_factorR == 0)
        for(uint i_r = 1; i_r < N_r; i_r++)
            bin_writer.write((char *)&listR[i_r], 8);
    if(log_factorPh == 0)
        for(uint i_ph = 1; i_ph < N_ph; i_ph++)
            bin_writer.write((char *)&listPh[i_ph], 8);
    if(log_factorTh == 0)
        for(uint i_th = 1; i_th < N_th; i_th++)
            bin_writer.write((char *)&listTh[i_th], 8);

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        cout << "-> Writing binary spherical grid file: " << float(100.0 * double(i_r) / double(N_r))
             << "      \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_th = 0; i_th < N_th; i_th++)
            {
                for(uint i = 0; i < data_offset; i++)
                {
                    double tmp_data = grid_cells[i_r][i_ph][i_th]->getData(i);
                    bin_writer.write((char *)&tmp_data, 8);
                }
            }
        }
    }

    for(uint i = 0; i < data_offset; i++)
    {
        double tmp_data = center_cell->getData(i);
        bin_writer.write((char *)&tmp_data, 8);
    }

    bin_writer.close();

    cout << CLR_LINE;
    cout << "- Writing spherical grid file   : done" << endl;

    return true;
}

bool CGridSpherical::createArtificialGrid(string path)
{
    resetGridValues();

    line_counter = 1;
    char_counter = 0;

    data_offset = 9;
    data_len = 0;
    max_data = 9;
    data_pos_gd_list.push_back(0);
    data_pos_dt_list.push_back(1);
    data_pos_tg = 2;
    data_pos_mx = 3;
    data_pos_my = 4;
    data_pos_mz = 5;
    data_pos_vx = 6;
    data_pos_vy = 7;
    data_pos_vz = 8;

    min_len = 1e30;
    max_len = 2 * Rmax;

    string filename = path;

    dataID = GRID_ID_SPH;
    Rmin = 1; //* con_pc;
    Rmax = 4; //*con_pc;
    N_r = 3;
    N_ph = 4;
    N_th = 4;
    log_factorR = 0;
    log_factorPh = 0;
    log_factorTh = 0;

    ofstream bin_writer(filename.c_str(), ios::out | ios::binary);

    if(bin_writer.fail())
    {
        cout << "\nERROR: Cannot write to:\n" << endl;
        cout << filename;
        return false;
    }

    bin_writer.write((char *)&dataID, 2);
    bin_writer.write((char *)&data_offset, 2);

    ushort tmp_ids;

    tmp_ids = GRIDgas_dens;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDdust_temp;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDgas_temp;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDmx;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDmy;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDmz;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDvx;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDvy;
    bin_writer.write((char *)&tmp_ids, 2);

    tmp_ids = GRIDvz;
    bin_writer.write((char *)&tmp_ids, 2);

    bin_writer.write((char *)&Rmin, 8);
    bin_writer.write((char *)&Rmax, 8);
    bin_writer.write((char *)&N_r, 2);
    bin_writer.write((char *)&N_ph, 2);
    bin_writer.write((char *)&N_th, 2);
    bin_writer.write((char *)&log_factorR, 8);
    bin_writer.write((char *)&log_factorPh, 8);
    bin_writer.write((char *)&log_factorTh, 8);

    for(uint i_r = 0; i_r < N_r; i_r++)
    {
        if(i_r % 50 == 0)
            cout << "-> Creating tree: " << 100.0 * float(i_r) / float(N_r) << " [%]           \r" << flush;

        for(uint i_ph = 0; i_ph < N_ph; i_ph++)
        {
            for(uint i_th = 0; i_th < N_th; i_th++)
            {
                max_cells++;

                double tmp_data = 1.0 / double(i_r * i_r + 1);

                bin_writer.write((char *)&tmp_data, 8);

                tmp_data = 1;
                bin_writer.write((char *)&tmp_data, 8);

                tmp_data = i_ph + 10;
                bin_writer.write((char *)&tmp_data, 8);

                for(uint i = 3; i < 9; i++)
                {
                    tmp_data = double(10);
                    bin_writer.write((char *)&tmp_data, 8);
                }
            }
        }
    }

    double tmp_data = 1e-5;
    bin_writer.write((char *)&tmp_data, 8);

    tmp_data = 1;
    bin_writer.write((char *)&tmp_data, 8);

    tmp_data = 2;
    bin_writer.write((char *)&tmp_data, 8);

    for(uint i = 3; i < 9; i++)
    {
        tmp_data = double(1);
        bin_writer.write((char *)&tmp_data, 8);
    }

    bin_writer.close();

    cout << "min: " << min_gas_dens << "  max_dens: " << max_gas_dens << endl;
    cout << "Creating artificial sphere                  : done" << endl;
    cout << "Max cells: " << max_cells << endl;

    return true;
}

void CGridSpherical::printParameters()
{
    if(max_cells == 0)
        cout << "\nERROR: No tree parameters available! " << endl;
    else
    {
        cout << CLR_LINE;
        cout << "Spherical grid parameters (ID: " << getDataID() << ", data len.: " << getDataOffset()
             << ", Nr: " << N_r << ", Nph: " << N_ph << ", Nth: " << N_th << ")" << endl;
        cout << SEP_LINE;

        cout << "- Number of spherical cells     : " << max_cells << endl;
        printPhysicalParameters();
        cout << SEP_LINE;
    }
}

bool CGridSpherical::positionPhotonInGrid(photon_package * pp)
{
    uint dirID = pp->getDirectionID();
    uint rID = MAX_UINT, thID = MAX_UINT, phID = MAX_UINT;
    if(dirID < 6 && pp->getPositionCell() != 0)
    {
        cell_sp * tmp_cell = (cell_sp *)pp->getPositionCell();
        bool skip = false;

        rID = tmp_cell->getRID();
        thID = tmp_cell->getThID();
        phID = tmp_cell->getPhID();

        // Update index of next cell
        switch(dirID)
        {
            case 0:
                rID--;
                break;

            case 1:
                rID++;
                break;

            case 2:
                thID--;
                break;

            case 3:
                thID++;
                break;

            case 4:
                if(phID == 0)
                    phID += N_ph;
                phID--;
                break;

            case 5:
                phID++;
                if(phID >= N_ph)
                    phID -= N_ph;
                break;

            default:
                return false;
                break;
        }
    }

    Vector3D sp_pos = pp->getPosition().getSphericalCoord();

    if(Rmin > sp_pos.R())
    {
        pp->setPositionCell(center_cell);
        return true;
    }

    uint i_r = 0, i_ph = 0, i_th = 0;

    i_r = CMathFunctions::biListIndexSearch(sp_pos.R(), listR, N_r + 1);
    if(i_r == MAX_UINT)
        return false;

    if(N_ph > 1)
    {
        i_ph = CMathFunctions::biListIndexSearch(sp_pos.Phi(), listPh, N_ph + 1);
        if(i_ph == MAX_UINT)
            return false;
    }

    i_th = CMathFunctions::biListIndexSearch(sp_pos.Z(), listTh, N_th + 1);
    if(i_th == MAX_UINT)
        return false;

    pp->setPositionCell(grid_cells[i_r][i_ph][i_th]);

    return true;
}

bool CGridSpherical::goToNextCellBorder(photon_package * pp)
{
    cell_sp * tmp_cell = (cell_sp *)pp->getPositionCell();
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();

    bool hit = false;
    double min_length = 1e300;
    double tmp_length[4];
    uint dirID = MAX_UINT;

    uint rID = tmp_cell->getRID();

    if(rID == MAX_UINT)
    {
        double B = 2 * p * d;
        double C = p.sq_length() - Rmin * Rmin;
        double dscr = B * B - 4 * C;

        if(dscr >= 0)
        {
            dscr = sqrt(dscr);
            tmp_length[0] = (-B + dscr) / 2;
            tmp_length[1] = (-B - dscr) / 2;
        }
        else
        {
            tmp_length[0] = 1e200;
            tmp_length[1] = 1e200;
        }

        for(uint i = 0; i < 2; i++)
        {
            if(tmp_length[i] >= 0 && tmp_length[i] < min_length)
            {
                min_length = tmp_length[i];
                hit = true;
                dirID = 1;
            }
        }
    }
    else
    {
        // --- Radial cell borders ---

        double r1 = listR[rID];
        double r2 = listR[rID + 1];

        double p_sq = p.sq_length();
        double B = 2 * p * d;
        double B_sq = pow(B, 2);

        double C1 = p_sq - r1 * r1;
        double C2 = p_sq - r2 * r2;

        double dscr1 = B_sq - 4 * C1;
        double dscr2 = B_sq - 4 * C2;

        if(dscr1 >= 0)
        {
            dscr1 = sqrt(dscr1);
            tmp_length[0] = (-B + dscr1) / 2;
            tmp_length[1] = (-B - dscr1) / 2;
        }
        else
        {
            tmp_length[0] = 1e200;
            tmp_length[1] = 1e200;
        }

        if(dscr2 >= 0)
        {
            dscr2 = sqrt(dscr2);
            tmp_length[2] = (-B + dscr2) / 2;
            tmp_length[3] = (-B - dscr2) / 2;
        }
        else
        {
            tmp_length[2] = 1e200;
            tmp_length[3] = 1e200;
        }

        for(uint i = 0; i < 4; i++)
        {
            if(tmp_length[i] >= 0 && tmp_length[i] < min_length)
            {
                min_length = tmp_length[i];
                hit = true;
                dirID = uint(i / 2.0);
            }
        }

        // --- Theta cell borders ---
        uint thID = tmp_cell->getThID();

        double th1 = listTh[thID];
        double cos_th1 = cos(th1);
        double cos_th1_sq = cos_th1 * cos_th1;

        double A1 = d.Z() * d.Z() - cos_th1_sq;
        double B1 = d.Z() * p.Z() * (1 - cos_th1_sq) - cos_th1_sq * (d.X() * p.X() + d.Y() * p.Y());
        double C3 = p.Z() * p.Z() * (1 - cos_th1_sq) - cos_th1_sq * (p.X() * p.X() + p.Y() * p.Y());

        double dscr3 = B1 * B1 - A1 * C3;

        if(dscr3 >= 0)
        {
            dscr3 = sqrt(dscr3);
            tmp_length[0] = (-B1 + dscr3) / A1;
            tmp_length[1] = (-B1 - dscr3) / A1;
        }
        else
        {
            tmp_length[0] = 1e200;
            tmp_length[1] = 1e200;
        }

        for(uint i = 0; i < 2; i++)
        {
            if(tmp_length[i] >= 0 && tmp_length[i] < min_length)
            {
                if((p.Z() + d.Z() * tmp_length[i]) * cos_th1 > 0)
                {
                    min_length = tmp_length[i];
                    hit = true;
                    dirID = 2;
                }
            }
        }

        double th2 = listTh[thID + 1];
        double cos_th2 = cos(th2);
        double cos_th2_sq = cos_th2 * cos_th2;

        double A2 = d.Z() * d.Z() - cos_th2_sq;
        double B2 = d.Z() * p.Z() * (1 - cos_th2_sq) - cos_th2_sq * (d.X() * p.X() + d.Y() * p.Y());
        double C4 = p.Z() * p.Z() * (1 - cos_th2_sq) - cos_th2_sq * (p.X() * p.X() + p.Y() * p.Y());

        double dscr4 = B2 * B2 - A2 * C4;

        if(dscr4 >= 0)
        {
            dscr4 = sqrt(dscr4);
            tmp_length[0] = (-B2 + dscr4) / A2;
            tmp_length[1] = (-B2 - dscr4) / A2;
        }
        else
        {
            tmp_length[0] = 1e200;
            tmp_length[1] = 1e200;
        }

        for(uint i = 0; i < 2; i++)
        {
            if(tmp_length[i] >= 0 && tmp_length[i] < min_length)
            {
                if((p.Z() + d.Z() * tmp_length[i]) * cos_th2 > 0)
                {
                    min_length = tmp_length[i];
                    hit = true;
                    dirID = 3;
                }
            }
        }

        // --- Phi cell borders ---
        if(N_ph > 1)
        {
            double r = sqrt(p.sq_length());
            double rho = sqrt(p.X() * p.X() + p.Y() * p.Y());

            double sin_th = rho / r;
            double cos_th = p.Z() / r;

            uint phID = tmp_cell->getPhID();
            double ph1 = listPh[phID];
            double ph2 = listPh[phID + 1];

            double sin_ph1 = sin(ph1);
            double sin_ph2 = sin(ph2);

            double cos_ph1 = cos(ph1);
            double cos_ph2 = cos(ph2);

            Vector3D v_n1 = -Vector3D(-sin_ph1, cos_ph1, 0);
            Vector3D v_a1 = r * Vector3D(sin_th * cos_ph1, sin_th * sin_ph1, cos_th);

            double den1 = v_n1 * d;
            if(den1 != 0)
            {
                double num = v_n1 * (p - v_a1);
                double length = -num / den1;

                if(length >= 0 && length < min_length)
                {
                    min_length = length;
                    hit = true;
                    dirID = 4;
                }
            }

            Vector3D v_n2 = Vector3D(-sin_ph2, cos_ph2, 0);
            Vector3D v_a2 = r * Vector3D(sin_th * cos_ph2, sin_th * sin_ph2, cos_th);

            double den2 = v_n2 * d;
            if(den2 != 0)
            {
                double num = v_n2 * (p - v_a2);
                double length = -num / den2;

                if(length >= 0 && length < min_length)
                {
                    min_length = length;
                    hit = true;
                    dirID = 5;
                }
            }
        }
    }

    if(!hit)
    {
        cout << "\nERROR: Wrong cell border!                                   " << endl;
        return false;
    }

    double path_length = min_length + 1e-3 * min_len;
    pp->setPosition(p + d * path_length);
    pp->setTmpPathLength(path_length);
    pp->setDirectionID(dirID);
    return true;
}

bool CGridSpherical::updateShortestDistance(photon_package * pp)
{
    Vector3D tmp_pos;
    double min_dist, tmp_dist[6];

    double loc_x_min, loc_x_max, loc_y_min, loc_y_max, loc_z_min, loc_z_max;
    bool found = false;

    cell_oc * tmp_cell_pos = (cell_oc *)pp->getPositionCell();

    tmp_pos = pp->getPosition();

    loc_x_min = tmp_cell_pos->getXmin();
    loc_y_min = tmp_cell_pos->getYmin();
    loc_z_min = tmp_cell_pos->getZmin();

    loc_x_max = tmp_cell_pos->getXmax();
    loc_y_max = tmp_cell_pos->getYmax();
    loc_z_max = tmp_cell_pos->getZmax();

    min_dist = 1E200;

    tmp_dist[0] = abs(loc_x_max - tmp_pos.X());
    tmp_dist[1] = abs(tmp_pos.X() - loc_x_min);

    tmp_dist[2] = abs(loc_y_max - tmp_pos.Y());
    tmp_dist[3] = abs(tmp_pos.Y() - loc_y_min);

    tmp_dist[4] = abs(loc_z_max - tmp_pos.Z());
    tmp_dist[5] = abs(tmp_pos.Z() - loc_z_min);

    for(int i = 0; i < 6; i++)
    {
        if(min_dist > tmp_dist[i])
        {
            min_dist = tmp_dist[i];
            found = true;
        }
    }

    pp->setShortestDistance(min_dist);
    return found;
}

bool CGridSpherical::findStartingPoint(photon_package * pp)
{
    Vector3D p = pp->getPosition();
    Vector3D d = pp->getDirection();

    if(isInside(p))
        return positionPhotonInGrid(pp);

    double tmp_length[2];
    double min_length = 1e200;
    bool hit = false;

    double B = 2 * p * d;
    double C = p.sq_length() - Rmax * Rmax;
    double dscr = B * B - 4 * C;

    if(dscr >= 0)
    {
        dscr = sqrt(dscr);
        tmp_length[0] = (-B + dscr) / 2;
        tmp_length[1] = (-B - dscr) / 2;
    }
    else
    {
        tmp_length[0] = 1e200;
        tmp_length[1] = 1e200;
    }

    for(uint i = 0; i < 2; i++)
    {
        if(tmp_length[i] >= 0 && tmp_length[i] < min_length)
        {
            min_length = tmp_length[i];
            hit = true;
        }
    }

    if(!hit)
        return false;

    double path_length = min_length + 1e-3 * min_len;
    pp->setPosition(p + d * path_length);
    pp->setDirectionID(MAX_UINT);
    return positionPhotonInGrid(pp);
}

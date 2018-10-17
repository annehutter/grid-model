#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>

int file_exist(char *filename)
{
	FILE *file;
	if((file = fopen (filename, "rt"))){
		fclose(file);
		return 1;
	}else{
		return 0;
	}
}

int main()
{
    char string[10];
    int int_value;
    
    FILE *file;
    
    printf("This program allows you to interactively create the needed input file to run the grid model.\nIt will ask you for the details of your simulation run. After typing the response, press ENTER key to continue\n\nKeep in mind that this program is not failsafe, so consider your inputs carefully before pressing enter! If you recognize a typo, you can still edit the resulting input file manually.\n\n"); 
    printf("Do you have all data files ready?\n");
    //here also if you press any other key will wait till pressing ENTER
    scanf("%s",string);
    if(strcmp(string,"No") == 0 || strcmp(string,"no") == 0)
    {
        printf("Get your data files ready! You will need: \n - density grid(s)\n - clumping factor grid(s)\n - list of sources (x,y,z,#Nion,ID,fesc) or grid(s) with #ionizing photons\n - recombination lookup table if you want to compute HI fractions within ionized regions (available on git repository)\n");
        exit(0);
    }
    for(int i=0; i<10; i++) string[i] = '\0';
    printf("\n\n");
    
    int            simulationType;
    int            calc_ion_history;
    int            num_snapshots;
    char           redshift_file[1000];
    double         redshift_prev_snap;
    double         redshift;
    double         evol_time;
    int            snapshot_start;

    //Cosmology
    double         h;
    double         omega_b;
    double         omega_m;
    double         omega_l;
    double         sigma8;
    double         Y;
    
    //BubbleModel
    double         lin_scales;
    double         inc_log_scales;
    double         max_scale;
    int            ion_sphere;
        
    //Input
    int            grid_size;
    double         box_size;
    
    int            input_doubleprecision;
    int            inputfiles_comoving;

    char           igm_density_file[1000];
    int            dens_in_overdensity;
    double         mean_density;
    int            default_mean_density;

    char           igm_clump_file[1000];
    
    char           sources_file[1000];
    char           nion_file[1000];
    
    //Photoionization
    int            use_web_model;
    int            photHI_model;
    int            calc_mfp;
    double         photHI_bg;
    char           photHI_bg_file[1000];
    double         mfp;
    double         source_slope_index = 5.;

    //Recombinations
    int            recomb;
    int            const_recomb;
    int            calc_recomb;
    double         dnrec_dt;
    char           recomb_table[1000];
    double         zmin, zmax, dz;
    double         fmin, fmax, df;
    double         dcellmin, dcellmax, ddcell;
    
    //Helium
    int            solve_He;
    char           sources_HeI_file[1000];
    char           nion_HeI_file[1000];
    char           sources_HeII_file[1000];
    char           nion_HeII_file[1000];
    
    double         dnrec_HeI_dt;
    double         dnrec_HeII_dt;
    
    //Output
    char           out_XHII_file[1000];
    int            write_photHI_file;
    char           out_photHI_file[1000];
    char           out_XHeII_file[1000];
    char           out_XHeIII_file[1000];
    
    //Restart
    int            write_restart_file;
    char           restart_file[1000];
    double         walltime;
    
    
    //-------------------------------------------------------------------------------
    // GET INFORMATION FROM USER
    //-------------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------------
    // Simulation type
    //-------------------------------------------------------------------------------
    printf("\nTYPE OF SIMULATION\n");
    
    printf("What type of simulation would you like to run?\n1) compute the ionization field from a single snapshot at a fixed redshift assuming the sources emitted for t Myrs (type '0')\n2) compute the evolution of ionization fields from a single snapshot across a given redshift span (type '1')\n3) compute the evolution of ionization fields from multiple snapshots across redshift (type '2')");
    scanf("%d", &simulationType);
    
    if(simulationType == 0)
    {
        num_snapshots = 0;
        for(int i=0; i<1000; i++) redshift_file[i] = '\0';
        strcat(redshift_file, "None");
        redshift_prev_snap = 0;
        snapshot_start = -1;

        printf("\nREDSHIFT & EVOLUTION TIME\n");

        printf("Which redshift are you considering?\n");
        scanf("%lf", &redshift);
        
        printf("What is the evolution time of the ionized regions you want to consider? Please specify in Myrs (double)\n");
        scanf("%lf", &evol_time);
        printf("evolutionTime = %e\n", evol_time);
        if(evol_time <= 0. || evol_time > 1.4e4)
        {
            printf("This time is not valid, it is smaller than 0. or larger than the age of the Universe\nPlease revise!\n");
            exit(0);
        }
    }
    else if(simulationType == 1)
    {
        for(int i=0; i<1000; i++) redshift_file[i] = '\0';
        strcat(redshift_file, "None");        
        evol_time = 0.;
        snapshot_start = -1;

        printf("\nREDSHIFT INTERVALS\n");

        printf("How many steps (including the initial AND final redshift, as well as all redshifts at which you provide input files) do you want to consider?\n");
        scanf("%d", &num_snapshots);
        
        printf("At which redshift starts your simulation? (double)\n");
        scanf("%lf", &redshift_prev_snap);
        
        printf("And at which redshift does it end? (double)\n");
        scanf("%lf", &redshift);
    }
    else if(simulationType == 2)
    {
        redshift = 0.;
        redshift_prev_snap = 0.;
        evol_time = 0.;
        
        printf("How many steps (including the initial AND final redshift, as well as all redshifts at which you provide input files) do you want to consider?\n");
        scanf("%d", &num_snapshots);
        
        printf("You need to provide a file which lists from high to low redshifts in each line:\n <redshift> <0/1>\nwhereas 0 means no input files at this redshift and 1 means input files are provided at this redshift. Please specify the address of this file now\n");
        scanf("%s", redshift_file);
        
        if(file_exist(redshift_file) == 0)
        {
            printf("WARNING:  Don't forget to create that file or check the address!\n");
        }
 
        printf("Please specify the first snapshot number of your input files (density, clumping factor, sources), i.e. _00i\n");
        scanf("%d", &snapshot_start);
        if(snapshot_start <= 0)
        {
            simulationType = 3;
            snapshot_start = -1;
        }
    }
    else
    {
        printf("This is not a valid input.\n");
        exit(EXIT_FAILURE);
    }
    
    
    //-------------------------------------------------------------------------------
    // Filtering characteristics
    //-------------------------------------------------------------------------------
    
    printf("\nFILTERING CHARACTERISTICS\n");

    printf("\nThis code uses various filter scales to compute the sizes of the ionized regions. Until which scale do you want to increase the filtering scale linearly? (in h^-1 Mpc) (double)\n");
    scanf("%lf", &lin_scales);
    printf("And how large should be the first increment on a log scale? (in h^-1 Mpc) (double)\n");
    scanf("%lf", &inc_log_scales);
    printf("And what should be the largest filtering scale? (in h^-1 Mpc) (double)\n");
    scanf("%lf", &max_scale);
    ion_sphere = 0;
    printf("And do you want to flag the entire sphere as ionized (1) or just the central cell (0)?\n");
    scanf("%d", &ion_sphere);
    
    //-------------------------------------------------------------------------------
    // OPTIONS
    //-------------------------------------------------------------------------------
    printf("\nOPTIONS\n");

    printf("\nThis code has multiple options on which properties to include or how to calculate those\n");
    
    
    //-------------------------------------------------------------------------------
    // Density
    //-------------------------------------------------------------------------------
    printf("\nDENSITY\n");

    printf("Do you want the default mean density, i.e. particle number density at the given redshift? Yes = 1, No = 0\n");
    scanf("%d", &default_mean_density);
    
    if(default_mean_density == 1)
    {
        mean_density = 0.;
    }
    else if(default_mean_density == 0)
    {
        printf("Please state the assumed mean density at z=0 (double):\n");
        scanf("%lf", &mean_density);
    }
    else
    {
        printf("This input is not valid\n");
        exit(0);
    }
    
    //-------------------------------------------------------------------------------
    // Helium
    //-------------------------------------------------------------------------------
    printf("\nHELIUM\n");

    printf("Do you want to consider also helium? Yes = 1, No = 0\n");
    scanf("%d", &solve_He);
    
    //-------------------------------------------------------------------------------
    // Recombinations
    //-------------------------------------------------------------------------------
    printf("\nRECOMBINATIONS\n");
    
    printf("Do you want to consider recombinations? Yes = 1, No = 0\n");
    scanf("%d", &recomb);
    
    if(recomb == 1){
        printf("Do you want to use the web model to compute the number of recombinations (1) or assume a constant number of recombinations (0)? \n");
        scanf("%d", &calc_recomb);
        
        for(int i=0; i<1000; i++) recomb_table[i] = '\0';
        if(calc_recomb == 1)
        {
            printf("Do you want to compute the number of recombinations approximately from the rate equations (1) or using the subgrid density model in Miralda_Escude 2000 (2)?\n");
            scanf("%d", &calc_recomb);
            if(calc_recomb == 2)
            {
                printf("Please type in the directory where the file with the nrec values is located. The file can be obtained from the github repository.\n");
                scanf("%s", recomb_table);
            }
            
            const_recomb = 0;
            dnrec_dt = 0.;
            dnrec_HeI_dt = 0.;
            dnrec_HeII_dt = 0.;
        }
        else
        {
            const_recomb = 1;
            
            printf("Which rate of recombinations for HII (dnrec/dt) do you want to assume (in Myrs^-1)?\n");
            scanf("%lf", &dnrec_dt);
            
            if(solve_He == 1)
            {
                printf("You included helium into your calculations:\n");
                printf("Which rate of recombinations for HeII (dnrec/dt) do you want to assume (in Myrs^-1)?\n");
                scanf("%lf", &dnrec_HeI_dt);
                printf("Which rate of recombinations for HeIII (dnrec/dt) do you want to assume (in Myrs^-1)?\n");
                scanf("%lf", &dnrec_HeII_dt);
            }
        }
    }
    else
    {
        const_recomb = 0;
        dnrec_dt = 0.;
        dnrec_HeI_dt = 0.;
        dnrec_HeII_dt = 0.;
    }
    
    //-------------------------------------------------------------------------------
    // Web model
    //-------------------------------------------------------------------------------
    printf("\nWEB MODEL (INCLUDING PHOTOIONIZATION)\n");
        
    printf("Do you want to use the web model (assume/compute a photoionization background & recombinations)? Yes = 1, No = 0\n");
    scanf("%d", &use_web_model);
    if(use_web_model == 0)
    {
        photHI_model = 0;
        calc_mfp = 0;
        calc_recomb = 0;
        
        for(int i=0; i<1000; i++) 
        {
            photHI_bg_file[i] = '\0';
            recomb_table[i] = '\0';
        }
        strcat(photHI_bg_file, "None");
    }
    else if(use_web_model == 1)
    {
        printf("\nPHOTOIONIZATION RATE\n");

        printf("Do you want to assume a constant photoionization field? Yes = 0, No = 1\n");
        scanf("%d",&photHI_model);
        if(photHI_model == 0)
        {
            photHI_model = 0;
            calc_mfp = 0;
            strcat(photHI_bg_file,"None");
            
            printf("Which value do you want to assume for the photoionization background (in s^-1)?\n");
            scanf("%lf", &photHI_bg);
        }
        else if(photHI_model > 0)
        {
            printf("How do you want to calculate the photoionization rate? Assuming a mean free path and apply a r^-2 kernel (1), or derive it from the number of ionizing photons in the ionized regions given the largest filtering scale (2)?\n");
            scanf("%d",&photHI_model);

            if(photHI_model == 1)
            {
                printf("Do you want to calculate the mean free path from the size of the ionized regions and Miralda-Escude 2000 (type 1) or set a value (type 0)?\n");
                scanf("%d", &calc_mfp);
                
                if(calc_mfp == 1)
                {
                    mfp = 0.;
                }
                else
                {
                    printf("Which mean free path in the ionized medium to you want to assume (in Mpc)?\n");
                    scanf("%lf", &mfp); 
                }
                
//                 printf("Please provide a list of the photoionization values at different redshifts (z,photIonHI,photHeatHI,Q). Specify the address of this file:\n");
//                 scanf("%s", photHI_bg_file);
                photHI_bg = 0.;
            }
            else if(photHI_model == 2)
            {
                mfp = 0.;
                strcat(photHI_bg_file, "None");
                photHI_bg = 0.;
            }
            else
            {
                printf("This input is not valid\n");
                exit(0);
            }
            
            printf("Which slope index do you want to assume for the spectrum of the ionizing sources: f ~ nu^-alpha. Type alpha\n");
            scanf("%lf", &source_slope_index);
        }
        else
        {
            printf("This input is not valid\n");
            exit(0);
        }
    }
    else
    {
        printf("This input is not valid\n");
        exit(0);
    }
    
    strcat(recomb_table, "nrec_values_batch_z3_30_0.01_f-9_9_0.1_d-4_4_0.1.dat");
    if(file_exist(recomb_table) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");
    zmin = 3.;
    zmax = 30.;
    dz = 0.01;
    fmin = -9.;
    fmax = 9.;
    df = 0.1;
    dcellmin = -4.;
    dcellmax = 4.;
    ddcell = 0.1;
    
    //-------------------------------------------------------------------------------
    // Cosmology
    //-------------------------------------------------------------------------------
    printf("\nCOSMOLOGY\n");

    printf("Which Cosmology do you want to assume?\n h = ");
    scanf("%lf", &h);
    printf(" omega_b = ");
    scanf("%lf", &omega_b);
    printf(" omega_m = ");
    scanf("%lf", &omega_m);
    printf(" omega_l = ");
    scanf("%lf", &omega_l);
    printf(" sigma8 = ");
    scanf("%lf", &sigma8);
    printf(" helium MASS fraction Y = ");
    scanf("%lf", &Y);
    
    //-------------------------------------------------------------------------------
    // Input parameters
    //-------------------------------------------------------------------------------
    printf("\nINPUT PARAMETERS\n");
    
    printf("What is the size of your input grids? \n");
    scanf("%d", &grid_size);
    printf("And what is the size of your simulation box in h^-1 Mpc?\n");
    scanf("%lf", &box_size);
    
    printf("Are your file in double (1) or single (0) precision?\n");
    scanf("%d", &input_doubleprecision);
    printf("Are your input grids comoving? Yes = 1, No = 0. Currently only comoving input files are supported -> Type '1' \n");
    scanf("%d", &inputfiles_comoving);
    if(inputfiles_comoving != 1)
    {
        printf("Sorry this does not work then.\n");
        exit(0);
    }
    
    //-------------------------------------------------------------------------------
    // Input grid files
    //-------------------------------------------------------------------------------
    printf("\nINPUT GRID FILES\n");
    
    printf("In the following type in the addresses of the input files.\n density file: ");
    scanf("%s", igm_density_file);
    if(file_exist(igm_density_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");

    printf(" Is the density given in overdensity (1) or in number densities (0)?\n");
    scanf("%d", &dens_in_overdensity);
    
    printf(" clumping factor file: ");
    scanf("%s", igm_clump_file);
    if(file_exist(igm_clump_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");

    printf("Are you providing a source list (0) or a grid with #ionizing photons (1) ?\n");
    scanf("%d", &int_value);
    if(int_value == 0)
    {
        printf("source list for HI ionizing photons: ");
        scanf("%s", sources_file);
        if(file_exist(sources_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");
        
        for(int i=0; i<1000; i++) nion_file[i] = '\0';
        strcat(nion_file, "None");
    }
    else
    {
        printf("Nion for HI grid: ");
        scanf("%s", nion_file);
        if(file_exist(nion_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");

        for(int i=0; i<1000; i++) sources_file[i] = '\0';
        strcat(sources_file, "None");
    }
    
    if(solve_He == 1)
    {
        if(int_value == 0)
        {
            printf("source list for HeI ionizing photons: ");
            scanf("%s", sources_HeI_file);
            if(file_exist(sources_HeI_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");
            
            for(int i=0; i<1000; i++) nion_HeI_file[i] = '\0';
            strcat(nion_HeI_file, "None");
            
            printf("source list for HeII ionizing photons: ");
            scanf("%s", sources_HeII_file);
            if(file_exist(sources_HeII_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");
            
            for(int i=0; i<1000; i++) nion_HeII_file[i] = '\0';
            strcat(nion_HeII_file, "None");
        }
        else
        {
            printf("Nion for HeI grid: ");
            scanf("%s", nion_HeI_file);
            if(file_exist(nion_HeI_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");
            
            for(int i=0; i<1000; i++) sources_HeI_file[i] = '\0';
            strcat(sources_HeI_file, "None");
            
            printf("Nion for HeII grid: ");
            scanf("%s", nion_HeII_file);
            if(file_exist(nion_HeI_file) == 0) printf("WARNING:  Don't forget to create that file or check the address!\n");
            
            for(int i=0; i<1000; i++) sources_HeII_file[i] = '\0';
            strcat(sources_HeII_file, "None");
        }
    }
    else
    {
        for(int i=0; i<1000; i++) nion_HeI_file[i] = '\0';
        for(int i=0; i<1000; i++) sources_HeI_file[i] = '\0';
        for(int i=0; i<1000; i++) nion_HeII_file[i] = '\0';
        for(int i=0; i<1000; i++) sources_HeII_file[i] = '\0';

        strcat(sources_HeI_file, "None");
        strcat(sources_HeII_file, "None");
        strcat(nion_HeI_file, "None");
        strcat(nion_HeII_file, "None");
    }
    
    //-------------------------------------------------------------------------------
    // Output files
    //-------------------------------------------------------------------------------
    printf("\nOUPUT FILES\n");
    
    printf("Finally the addresses of the output files: \n");
    printf("XHII grid: ");
    scanf("%s", out_XHII_file);
    
    printf("Do you want to output the HI photoionization fields? Yes = 1, No = 0\n");
    scanf("%d", &write_photHI_file);
    if(write_photHI_file == 1)
    {
        printf("HI Photoionization grid: ");
        scanf("%s", out_photHI_file);
    }
    else
    {
        for(int i=0; i<1000; i++) out_photHI_file[i] = '\0';
        strcat(out_photHI_file, "None");
    }

    if(solve_He == 1)
    {
        printf("XHeII grid:");
        scanf("%s", out_XHeII_file);
        printf("XHeIII grid:");
        scanf("%s", out_XHeIII_file);
    }
    else
    {
        for(int i=0; i<1000; i++) out_XHeII_file[i] = '\0';
        strcat(out_XHeII_file, "None");
        for(int i=0; i<1000; i++) out_XHeIII_file[i] = '\0';
        strcat(out_XHeIII_file, "None");
    }
    
    //-------------------------------------------------------------------------------
    // Restart option
    //-------------------------------------------------------------------------------
    
    printf("\nRESTART OPTIONS\n");
    
    printf("Do you want to store restart files, from which the simulation could be resumed? Yes=1, No=0\n");
    scanf("%d", &write_restart_file);
    if(write_restart_file == 1)
    {
        printf("basename for restart files:");
        scanf("%s", restart_file);
        
        printf("If the simulation should be only run for a certain time and restart files written before, provide the run time in minutes, otherwise tupe '0' and restart files are written after each output.\n");
        scanf("%lf", &walltime);
    }
    else
    {
        for(int i=0; i<1000; i++) restart_file[i] = '\0';
        strcat(restart_file, "None");
        walltime = 0.;
    }
    
    //-------------------------------------------------------------------------------
    // WRITE FILE
    //-------------------------------------------------------------------------------
    
    printf("Please state the address and filename of the inputfile you want to generate:\n");
    scanf("%s", string);
        
    file = fopen(string, "wt");
    
    fprintf(file, "[General]\n");
    fprintf(file, "calcIonHistory = %d\n", calc_ion_history);
    fprintf(file, "numSnapshots = %d\n", num_snapshots);
    fprintf(file, "redshiftFile = %s\n", redshift_file);
    fprintf(file, "redshift_prevSnapshot = %f\n", redshift_prev_snap);
    fprintf(file, "finalRedshift = %f\n", redshift);
    fprintf(file, "evolutionTime = %f\n\n", evol_time);
    
    fprintf(file, "size_linear_scale = %f\n", lin_scales);
    fprintf(file, "first_increment_in_logscale = %f\n", inc_log_scales);
    fprintf(file, "max_scale = %f\n", max_scale);
    fprintf(file, "useIonizeSphereModel = %d\n\n", ion_sphere);
    
    fprintf(file, "useDefaultMeanDensity = %d\n\n", default_mean_density);
    
    fprintf(file, "useWebModel = %d\n", use_web_model);
    fprintf(file, "photHImodel = %d\n", photHI_model);
    fprintf(file, "calcMeanFreePath = %d\n", calc_mfp);
    fprintf(file, "constantRecombinations = %d\n", const_recomb);
    fprintf(file, "calcRecombinations = %d\n\n", calc_recomb);
    
    fprintf(file, "solveForHelium = %d\n\n\n", solve_He);
    
    
    fprintf(file, "[Input]\n");
    fprintf(file, "gridsize = %d\n", grid_size);
    fprintf(file, "boxsize = %f\n\n", box_size);

    fprintf(file, "inputFilesAreInDoublePrecision = %d\n", input_doubleprecision);
    fprintf(file, "inputFilesAreComoving = %d\n\n", inputfiles_comoving);
    
    fprintf(file, "inputIgmDensityFile = %s\n", igm_density_file);
    fprintf(file, "densityInOverdensity = %d\n", dens_in_overdensity);
    fprintf(file, "meanDensity = %e\n\n", mean_density);
    
    fprintf(file, "inputIgmClumpFile = %s\n", igm_clump_file);
    
    fprintf(file, "inputSourcesFile = %s\n", sources_file);
    fprintf(file, "inputNionFile = %s\n\n\n", nion_file);

    
    fprintf(file, "[Output]\n");
    fprintf(file, "output_XHII_file = %s\n", out_XHII_file);
    fprintf(file, "write_photHI_file = %d\n", write_photHI_file);
    fprintf(file, "output_photHI_file = %s\n\n\n", out_photHI_file);
    
    
    fprintf(file, "[Cosmology]\n");
    fprintf(file, "hubble_h = %f\n", h);
    fprintf(file, "omega_b = %f\n", omega_b);
    fprintf(file, "omega_m = %f\n", omega_m);
    fprintf(file, "omega_l = %f\n", omega_l);
    fprintf(file, "sigma8 = %f\n", sigma8);
    fprintf(file, "Y = %f\n\n\n", Y);
    
    
    fprintf(file, "[Photoionization]\n");
    fprintf(file, "photHI_bg_file = %s\n", photHI_bg_file);
    fprintf(file, "photHI_bg = %e\n", photHI_bg);
    fprintf(file, "meanFreePathInIonizedMedium = %e\n", mfp);
    fprintf(file, "sourceSlopeIndex = %e\n\n\n", source_slope_index);
    
    
    fprintf(file, "[Recombinations]\n");
    fprintf(file, "dnrec_dt = %e\n", dnrec_dt);
    fprintf(file, "recombinationTable = %s\n", recomb_table);
    fprintf(file, "zmin = %f\n", zmin);
    fprintf(file, "zmax = %f\n", zmax);
    fprintf(file, "dz = %f\n", dz);
    fprintf(file, "fmin = %f\n", fmin);
    fprintf(file, "fmax = %f\n", fmax);
    fprintf(file, "df = %f\n", df);
    fprintf(file, "dcellmin = %f\n", dcellmin);
    fprintf(file, "dcellmax = %f\n", dcellmax);
    fprintf(file, "ddcell = %f\n\n", ddcell);
    
    
    fprintf(file, "[Helium]\n");
    fprintf(file, "inputSourcesHeIFile = %s\n", sources_HeI_file);
    fprintf(file, "inputNionHeIFile = %s\n", nion_HeI_file);
    fprintf(file, "inputSourcesHeIIFile = %s\n", sources_HeII_file);
    fprintf(file, "inputNionHeIIFile = %s\n\n", nion_HeII_file);
    
    fprintf(file, "dnrec_HeI_dt = %e\n", dnrec_HeI_dt);
    fprintf(file, "dnrec_HeII_dt = %e\n\n", dnrec_HeII_dt);

    fprintf(file, "output_XHeII_file = %s\n", out_XHeII_file);
    fprintf(file, "output_XHeIII_file = %s\n", out_XHeIII_file);

    fclose(file);
    
    printf("\n\nYour inputfile has been written. Now you can use it to run the gridmodel with:\n./cifog <INPUTFILE>\n");
    
    char s[] = { 0xf0, 0x9f, 0x98, 0x8e, 0 };

    printf("\n         %s\n\n", s);
    
    return 0;
}

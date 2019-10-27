//Author: Aliaksei Chareshneu
//Year: 2018
//Purpose: Program should draw a graph of accessability of residues in protein structure and simultaneously the graph will show hydrophobicity of each aminoacid. Additional purpose was to implement some kind of scaling algorithm to support several resolutions.
/*
Functions: 
1. Reads a configuration file, which name is specified via the command line, and stores information regarding the PDB file to be processed and regarding the window size;
2. Reads the requested PDB file (only standart residues from ATOM records);
3. Finds C-alpha atoms of each aminoacid;
4. Calculates the residue-residue distances as distances between C-alpha atoms;
5. Calculates value of accessability of each residue as a number of residues which are closer than 14 Angstroms from the given residue;
5. Displays a curve of accessability;
6. Displays a histogram of hydrophobicity;
7. Displays a color coded marks at the bottom part of the graph together with one-letter codes of aminoacids and sequence numbers;
8. If the structure will have more than 50 residues, the graph will be divided to several separate graphs, each displaying maximum 50 residues;
9. Displays legend of color codes, title and the name of PDB file read;

Notification: Scaling has been implemented, but, to be honest, it is not perfect, so it would be better to test the program whithin the resolution range of [800; 1100] for both width and height (e.g. 800 x 1100, 1100 x 800, 800 x 800, 1100 x 1100, or any other value between 800 and 1100). The reason for such constraints is that currently the program in its scaling part partially rests on empirical values, and only partially on some variables, derived from width and height values, extracted from a configuration file. I have an idea of some clever scaling function or algorithm here, but it would require almost complete rewriting of significant part of graph drawing function.

Command line parameters:
The user is supposed to supply only one argument via the command line - name of the configuration file to be read. Otherwise, the program will report an error and ask the user to supply the name of the configuration file via the command line.

Format of the configuration file:

INPUT_FILE = PDB_filename
WINDOW_SIZE = width, height

There could be an arbitrary number of spaces before or after "=" sign.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <g2.h>
#include <g2_X11.h>
#include <math.h>
#include <errno.h>

/*maximum array size definition*/
#define MAX_ATOMS 100000
#define BUF_SIZE 500
#define MAXSIZE 10000
#define DEFAULT_OUTPUT "/home/aaxx/Desktop/2018_PhD/Study/C_language/c_test/output_final_1.pdb"
#define RESD_TYPES_COUNT 21

/*structures definition and variables initialization*/

    typedef struct
    {
        char record_name[7];
        int atom_number;
        char atom_name[5];
        char alt_location_indicator;
        char residue_name[4];
        char chain_id;
        int residue_seq_number;
        char insertion_code;
        double x;
        double y;
        double z;
        double occupancy;
        double temperature_factor;
        char element_symbol[3];
        char formal_charge[3];
        //when an input file has missed records (i.e. atom numbers increased by value, greater than 1), it is necessary to have such a correction. More information below
        int correction;
    } ATOM;
    
    typedef struct
    {
        int first_atom;
        int last_atom;
        int residue_seq_number_r;
        char residue_name_r[4];
        int atom_c_alpha;
        //The following is necessary to figure out the real position of first and last atom of each residue in the array of atoms (i.e., number of record, starting from 0, in the original PDB file) in case it has missed records
        int first_atom_correction;
        int last_atom_correction;
        //accessability here is a number of residues which are less than 14 Angstroms from the given residues/closer than 14 Angstroms to the given residue
        int accessability;
        int residue_type;
    } RESIDUE;
    
    typedef struct
    {
        char code3[4]; //three-letter code
        char code1;    //one-letter code
        double color_r;
        double color_g;
        double color_b;
        double hydrophobicity;
    } RESIDUE_TYPE;

/*initialization of global variables*/
    ATOM atoms[MAX_ATOMS];
    RESIDUE residues[MAX_ATOMS];
    
    RESIDUE_TYPE residue_types[RESD_TYPES_COUNT] = {{"UNK", 'X', 153/255.0, 153/255.0, 153/255.0, 0.0}, {"ALA", 'A', 204/255.0, 255/255.0, 255/255.0, 1.8}, {"ARG", 'R', 230/255.0, 6/255.0, 6/255.0, -4.5}, {"ASN", 'N', 255/255.0, 153/255.0, 0/255.0, -3.5}, {"ASP", 'D', 255/255.0, 204/255.0, 153/255.0, -3.5}, {"CYS", 'C', 0/255.0, 255/255.0, 255/255.0, 2.5}, {"GLN", 'Q', 255/255.0, 102/255.0, 0/255.0, -3.5}, {"GLU", 'E', 255/255.0, 204/255.0, 0/255.0, -3.5}, {"GLY", 'G', 0/255.0, 255/255.0, 0/255.0, -0.4}, {"HIS", 'H', 255/255.0, 255/255.0, 153/255.0, -3.2}, {"ILE", 'I', 0/255.0, 0/255.0, 128/255.0, 4.5}, {"LEU", 'L', 51/255.0, 102/255.0, 255/255.0, 3.8}, {"LYS", 'K', 198/255.0, 6/255.0, 0/255.0, -3.9}, {"MET", 'M', 153/255.0, 204/255.0, 255/255.0, 1.9}, {"PHE", 'F', 0/255.0, 204/255.0, 255/255.0, 2.8}, {"PRO", 'P', 255/255.0, 255/255.0, 0/255.0, -1.6}, {"SER", 'S', 204/255.0, 255/255.0, 153/255.0, -0.8}, {"THR", 'T', 0/255.0, 255/255.0, 153/255.0, -0.7}, {"TRP", 'W', 204/255.0, 153/255.0, 255/255.0, -0.9}, {"TYR", 'Y', 204/255.0, 255/255.0, 204/255.0, -1.3}, {"VAL", 'V', 0/255.0, 0/255.0, 255/255.0, 4.2}};
    
    int atom_count = 0;
    int atom_count_scout = 0;
    int residue_count = 0;
    
    char buf[BUF_SIZE] = "";
    char s[30] = "";
    
    char input_to_reading[BUF_SIZE] = "";
    //char input_to_writing[BUF_SIZE] = "";
    int width_from_configuration_file = 0;
    int height_from_configuration_file = 0;
    double distance_from_configuration_file = 0;
    //The following is an array for residue numbers
    int numbers_of_residues_to_be_read[MAXSIZE] = {0};
    //The following is a total number of residues which are to be read
    int total_number_of_residues_to_be_read = 0;
    //The following is a total number of residues in the input PDB file, which name is specified in the configuration file
    double N_residues = 1;
    //It is necessary to initialize dev2 variable as a global variable, since it will appear in main function. The reason for such decision is that there are two ways to enter the configuration file name - as a command line argument (by default) and, in case it was not supplied as a command line argument, the program will ask to supply it separately. So in that case it is necessary to keep the opened window with graph and histogram even when the user enters something (configuration file name) separately and hits Enter. That is why dev2 is declared as a global variable. More information is in the end of the program in main section and after it.
    int dev2 = 0;

// opening the file for reading    
int file_opening_and_reading(char input[MAXSIZE])
{
    FILE *f1 = NULL;
    /*file opening for reading*/    
    f1 = fopen(input, "r");

    /*checking whether the opening was successful*/    
    if (f1 == NULL)
    {
        printf("Cannot open the input file!\n");
        if (errno != 0)
            printf("Error explanation: %s\n", strerror(errno));
        return 1;
    }
    
    /*reading the file and checking for array overflow*/
    while (feof(f1) == 0)
        {
            memset(buf, '\0', BUF_SIZE);
            if (fgets(buf, BUF_SIZE, f1) == NULL)
                break;
            if (atom_count >= MAX_ATOMS)
            {
                printf("Size of the array is too small!11\n");
                break;
            }
            if (strncmp(buf, "ATOM", 4) == 0)
            {
                //record name reading
                strncpy(atoms[atom_count].record_name, buf, 6);
                atoms[atom_count].record_name[6] = '\0';
                //atom number reading
                strncpy(s, buf + 6, 5);
                s[5] = '\0';
                sscanf(s, "%d", &atoms[atom_count].atom_number);
                //atom name reading
                strncpy(atoms[atom_count].atom_name, buf + 12, 4);
                atoms[atom_count].atom_name[4] = '\0';
                //alternative location reading
                atoms[atom_count].alt_location_indicator = buf[16];
                //residue name reading
                strncpy(atoms[atom_count].residue_name, buf + 17, 3);
                atoms[atom_count].residue_name[3] = '\0';
                //chain ID reading
                atoms[atom_count].chain_id = buf[21];
                //residue sequence number reading
                strncpy(s, buf + 22, 4);
                s[4] = '\0';
                sscanf(s, "%d", &atoms[atom_count].residue_seq_number);
                //code for insertion of residues reading
                atoms[atom_count].insertion_code = buf[26];
                //x coordinate reading
                strncpy(s, buf + 30, 8);
                s[8] = '\0';
                sscanf(s, "%lf", &atoms[atom_count].x);
                //y coordinate reading
                strncpy(s, buf + 38, 8);
                s[8] = '\0';
                sscanf(s, "%lf", &atoms[atom_count].y);
                //z coordinate reading
                strncpy(s, buf + 46, 8);
                s[8] = '\0';
                sscanf(s, "%lf", &atoms[atom_count].z);
                //occupancy reading
                strncpy(s, buf + 54, 6);
                s[6] = '\0';
                sscanf(s, "%lf", &atoms[atom_count].occupancy);
                //temperature factor reading
                strncpy(s, buf + 60, 6);
                s[6] = '\0';
                sscanf(s, "%lf", &atoms[atom_count].temperature_factor);
                //element symbol reading
                strncpy(atoms[atom_count].element_symbol, buf + 76, 2);
                atoms[atom_count].element_symbol[2] = '\0';
                //formal charge on atoms reading
                strncpy(atoms[atom_count].formal_charge, buf + 78, 2);
                atoms[atom_count].formal_charge[2] = '\0';
                       
                atoms[atom_count].correction = atoms[atom_count].atom_number - atom_count;
                atom_count++;

            }
        }
        
        /*closing the file*/
        if (fclose(f1) == EOF)
        {
            printf("Error closing the file %s\n", input);
            if (errno != 0)
                printf("Explanation: %s\n", strerror(errno));
            return 1;
        }
        f1 = NULL;
        return 0;
}

/*opening the file for writing*/
int file_opening_and_writing(char output[MAXSIZE])
{
    int i = 0;
    int j = 0;
    
    FILE *f1w = NULL;
    
    f1w = fopen(output, "w");
    int atom_count2 = 0;
    
    if (f1w == NULL)
    {
        printf("Cannot open the output file!\n");
        if (errno != 0)
            printf("Error explanation: %s\n", strerror(errno));
        return 1;
    }
    for (i = 0; i < total_number_of_residues_to_be_read; i++)
    {
        j = numbers_of_residues_to_be_read[i] - 1;        
        for (atom_count2 = (residues[j].first_atom - residues[j].first_atom_correction); atom_count2 <= (residues[j].last_atom - residues[j].last_atom_correction); atom_count2++)
        {
            if (atom_count2 >= MAX_ATOMS)
            {
                printf("Size of the array is too small!11\n");
                break;
            }

            fprintf(f1w, "%-6.6s%5d %-4.4s%c%-3.3s %c%4d%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2.2s%2.2s\n", atoms[atom_count2].record_name, atoms[atom_count2].atom_number, atoms[atom_count2].atom_name, atoms[atom_count2].alt_location_indicator, atoms[atom_count2].residue_name, atoms[atom_count2].chain_id, atoms[atom_count2].residue_seq_number, atoms[atom_count2].insertion_code, atoms[atom_count2].x, atoms[atom_count2].y, atoms[atom_count2].z, atoms[atom_count2].occupancy, atoms[atom_count2].temperature_factor, atoms[atom_count2].element_symbol, atoms[atom_count2].formal_charge);
        }
    }
    if (fclose(f1w) == EOF)
        {
        printf("Error closing the file %s\n", output);
        if (errno != 0)
            printf("Explanation: %s\n", strerror(errno));
        return 1;
        }
    return 0;
}

void first_atom_scout()
{
    char str[30] = "";
    
    printf("\n");
    printf("====================================================\n");
    printf("List of first and last atom indexes for each residue\n");
    printf("====================================================\n");
    
    for (atom_count_scout = 0; atom_count_scout <= atom_count; atom_count_scout++)
    {
        if (atom_count_scout >= MAX_ATOMS)
        {
            printf("Size of the array is too small!11\n");
            break;
        }
        if (atoms[atom_count_scout].residue_seq_number > atoms[atom_count_scout - 1].residue_seq_number)
        {
//Adding new residues into the array of structures RESIDUE
            snprintf(str, MAXSIZE, "%d", atoms[atom_count_scout].atom_number);
            sscanf(str, "%d", &residues[residue_count].first_atom);
            snprintf(str, MAXSIZE, "%d", atoms[atom_count_scout-1].atom_number);
            sscanf(str, "%d", &residues[residue_count-1].last_atom);
            snprintf(str, MAXSIZE, "%d", atoms[atom_count_scout].residue_seq_number);
            sscanf(str, "%d", &residues[residue_count].residue_seq_number_r);
            snprintf(str, MAXSIZE, "%s", atoms[atom_count_scout].residue_name);
            sscanf(str, "%s", residues[residue_count].residue_name_r);
            //correction values for the first and the last atoms in case of missed records in the original PDB file
            snprintf(str, MAXSIZE, "%d", atoms[atom_count_scout].correction);
            sscanf(str, "%d", &residues[residue_count].first_atom_correction);
            snprintf(str, MAXSIZE, "%d", atoms[atom_count_scout-1].correction);
            sscanf(str, "%d", &residues[residue_count-1].last_atom_correction);
            
            //Adding values into residue_type variable in residues structure
            if (strcmp(residues[residue_count].residue_name_r, "UNK") == 0)
                residues[residue_count].residue_type = 0;
            else if (strcmp(residues[residue_count].residue_name_r, "ALA") == 0)
                residues[residue_count].residue_type = 1;
            else if (strcmp(residues[residue_count].residue_name_r, "ARG") == 0)
                residues[residue_count].residue_type = 2;
            else if (strcmp(residues[residue_count].residue_name_r, "ASN") == 0)
                residues[residue_count].residue_type = 3;
            else if (strcmp(residues[residue_count].residue_name_r, "ASP") == 0)
                residues[residue_count].residue_type = 4;
            else if (strcmp(residues[residue_count].residue_name_r, "CYS") == 0)
                residues[residue_count].residue_type = 5;
            else if (strcmp(residues[residue_count].residue_name_r, "GLN") == 0)
                residues[residue_count].residue_type = 6;
            else if (strcmp(residues[residue_count].residue_name_r, "GLU") == 0)
                residues[residue_count].residue_type = 7;
            else if (strcmp(residues[residue_count].residue_name_r, "GLY") == 0)
                residues[residue_count].residue_type = 8;
            else if (strcmp(residues[residue_count].residue_name_r, "HIS") == 0)
                residues[residue_count].residue_type = 9;
            else if (strcmp(residues[residue_count].residue_name_r, "ILE") == 0)
                residues[residue_count].residue_type = 10;
            else if (strcmp(residues[residue_count].residue_name_r, "LEU") == 0)
                residues[residue_count].residue_type = 11;
            else if (strcmp(residues[residue_count].residue_name_r, "LYS") == 0)
                residues[residue_count].residue_type = 12;
            else if (strcmp(residues[residue_count].residue_name_r, "MET") == 0)
                residues[residue_count].residue_type = 13;
            else if (strcmp(residues[residue_count].residue_name_r, "PHE") == 0)
                residues[residue_count].residue_type = 14;
            else if (strcmp(residues[residue_count].residue_name_r, "PRO") == 0)
                residues[residue_count].residue_type = 15;
            else if (strcmp(residues[residue_count].residue_name_r, "SER") == 0)
                residues[residue_count].residue_type = 16;
            else if (strcmp(residues[residue_count].residue_name_r, "THR") == 0)
                residues[residue_count].residue_type = 17;
            else if (strcmp(residues[residue_count].residue_name_r, "TRP") == 0)
                residues[residue_count].residue_type = 18;
            else if (strcmp(residues[residue_count].residue_name_r, "TYR") == 0)
                residues[residue_count].residue_type = 19;
            else if (strcmp(residues[residue_count].residue_name_r, "VAL") == 0)
                residues[residue_count].residue_type = 20;
            
            residue_count++;
            N_residues++;
        }
        
//dealing with the very last atom
        if ((residues[residue_count].first_atom == 0) && (residues[residue_count].last_atom == 0))
        {
            snprintf(str, MAXSIZE, "%d", atoms[atom_count-1].atom_number);
            sscanf(str, "%d", &residues[residue_count-1].last_atom);
        }
    }
    
    for (atom_count_scout = 0; atom_count_scout < residue_count; atom_count_scout++)
    {
        printf("Residue: %d %s, indexes of atoms: %d, %d\n", residues[atom_count_scout].residue_seq_number_r, residues[atom_count_scout].residue_name_r, residues[atom_count_scout].first_atom, residues[atom_count_scout].last_atom);
    }
}

void barycentre_calc()
{
    int i = 0;
    int j = 0;
    double x_sum = 0;
    double y_sum = 0;
    double z_sum = 0;
    double x_bary = 0;
    double y_bary = 0;
    double z_bary = 0;
    int H_count = 0;
    
    printf("\n");
    printf("====================================\n");
    printf("List of barycentres for each residue\n");
    printf("====================================\n");
    for (i = 0; i < (atom_count_scout); i++)
    {
        for (j = (residues[i].first_atom - residues[i].first_atom_correction); j <= (residues[i].last_atom - residues[i].last_atom_correction); j++)
        {
            
            if (strcmp(atoms[j].element_symbol, " H") != 0)
            {
                x_sum = x_sum + atoms[j].x;
                y_sum = y_sum + atoms[j].y;
                z_sum = z_sum + atoms[j].z;
            }
            else
                H_count++;
        }
        x_bary = (x_sum)/((residues[i].last_atom - residues[i].last_atom_correction) - (residues[i].first_atom - residues[i].first_atom_correction) - H_count + 1);
        y_bary = (y_sum)/((residues[i].last_atom - residues[i].last_atom_correction) - (residues[i].first_atom - residues[i].first_atom_correction) - H_count + 1);
        z_bary = (z_sum)/((residues[i].last_atom - residues[i].last_atom_correction) - (residues[i].first_atom - residues[i].first_atom_correction) - H_count + 1);
        printf("Residue %d %s, barycentre: %lf, %lf, %lf\n", residues[i].residue_seq_number_r, residues[i].residue_name_r, x_bary, y_bary, z_bary);
        x_sum = 0;
        y_sum = 0;
        z_sum = 0;
        H_count = 0;
    }  
}

void backbone_atoms()
{
    int i = 0;
    int j = 0;
    int k = 0;
    int q = 0;
    int str_backbone_numbers[BUF_SIZE] = {0};
    
    printf("\n");
    printf("==============================================\n");
    printf("List of backbone atom numbers for each residue\n");
    printf("==============================================\n");
    for (i = 0; i < (atom_count_scout); i++)
    {
        for (j = (residues[i].first_atom - residues[i].first_atom_correction); j <= (residues[i].last_atom - residues[i].last_atom_correction); j++)
        {
            if (strcmp(atoms[j].atom_name, " N  ") == 0 || strcmp(atoms[j].atom_name, " CA ") == 0 || strcmp(atoms[j].atom_name, " C  ") == 0 || strcmp(atoms[j].atom_name, " O  ") == 0)
            {
                str_backbone_numbers[k] = atoms[j].atom_number;
                k++;
            }
        }
        str_backbone_numbers[k] = '\0';
        printf("Residue %d %s, backbone atom numbers: ", residues[i].residue_seq_number_r, residues[i].residue_name_r);
        for (q = 0; q < k; q++)
        {
            if (q >= MAX_ATOMS)
            {
                printf("Size of the array is too small!11\n");
                break;
            }
            printf("%d ", str_backbone_numbers[q]);
        }
        printf("\n");
        memset(&str_backbone_numbers[0], 0, sizeof(str_backbone_numbers));
        k = 0;
    }
}

void C_alpha_scout()
{
    int i = 0;
    int j = 0;
    
    for (i = 0; i <= (atom_count_scout); i++)
    {
        for (j = (residues[i].first_atom - residues[i].first_atom_correction); j <= (residues[i].last_atom - residues[i].last_atom_correction); j++)
        {
            if (strcmp(atoms[j].atom_name, " CA ") == 0)
            {
                residues[i].atom_c_alpha = j;
                break;
            }
        }
    }
}

double C_alpha_distance_calc(int index_1, int index_2)
{
    double distance = 0.0;
    distance = sqrt((atoms[index_1].x - atoms[index_2].x)*(atoms[index_1].x - atoms[index_2].x) + (atoms[index_1].y - atoms[index_2].y)*(atoms[index_1].y - atoms[index_2].y) + (atoms[index_1].z - atoms[index_2].z)*(atoms[index_1].z - atoms[index_2].z));
    return distance;
}

void C_alpha_residues_comparing()
{
    int i = 0;
    int j = 0;
    int k = 0;
    
    printf("\n");
    printf("==============================================\n");
    printf("List of distances between C-alpha atoms of residues, which are less than 14 (in Angstroms)\n");
    printf("==============================================\n");
    
    for (i = 0; i < (atom_count_scout); i++)
    {
        for (j = 0; j < (atom_count_scout); j++)
        {
            if (j >= MAX_ATOMS)
            {
                printf("Size of the array is too small!11\n");
                break;
            }
            if (((C_alpha_distance_calc(residues[i].atom_c_alpha, residues[j].atom_c_alpha)) < 14) && ((C_alpha_distance_calc(residues[i].atom_c_alpha, residues[j].atom_c_alpha)) != 0) && (strncmp(atoms[(residues[i].first_atom - residues[i].first_atom_correction)].record_name, "ATOM", 4) == 0) && (strncmp(atoms[(residues[j].first_atom - residues[j].first_atom_correction)].record_name, "ATOM", 4) == 0))
            {
                //printf("The distance between %d %s residue and %d %s residue is %lf\n", residues[i].residue_seq_number_r, residues[i].residue_name_r, residues[j].residue_seq_number_r, residues[j].residue_name_r, C_alpha_distance_calc(residues[i].atom_c_alpha, residues[j].atom_c_alpha));
                residues[i].accessability = residues[i].accessability + 1;
                //residues[j].accessability = residues[j].accessability + 1;
                k++;
            }
        }
        printf("Accessability for %d %s residue is %d\n", residues[i].residue_seq_number_r, residues[i].residue_name_r, residues[i].accessability);
    }
}

int configuration_file_reading (char input[MAXSIZE])
{
    int i = 0;
    int delta = 0;
    char RESIDUE_LIST[BUF_SIZE] = "";
    int N_of_characters_read = 0;
    FILE *fc = NULL;
    /*file opening for reading*/    
    fc = fopen(input, "r");

    /*checking whether the opening was successful*/    
    if (fc == NULL)
    {
        printf("Cannot open the configuration file!\n");
        if (errno != 0)
            printf("Error explanation: %s\n", strerror(errno));
        return 1;
    }
    
    /*reading the file and checking for array overflow*/
    while (feof(fc) == 0)
        {
            memset(buf, '\0', BUF_SIZE);
            if (fgets(buf, BUF_SIZE, fc) == NULL)
                break;
            if (strncmp(buf, "INPUT_FILE", 10) == 0)
            {
                //input filename reading
                sscanf(buf, " INPUT_FILE = %s ", input_to_reading);
            }
            else if (strncmp(buf, "WINDOW_SIZE", 11) == 0)
            {
                //window size reading
                //it is not that usefull here, in fact, it is not used at all, but it may be usefull later
                sscanf(buf, " WINDOW_SIZE = %d , %d ", &width_from_configuration_file, &height_from_configuration_file);
            }
            else if (strncmp(buf, "DISTANCE", 8) == 0)
            {
                //distance reading
                sscanf(buf, " DISTANCE = %lf ", &distance_from_configuration_file);
            }
            else if (strncmp(buf, "RESIDUE_LIST", 12) == 0)
            {
                sscanf(buf, "%[RESIDUE_LIST = ] %n", RESIDUE_LIST, &N_of_characters_read);
                delta = N_of_characters_read;
                for (i = 0; i < MAXSIZE; i++)
                {
                    sscanf(buf + delta, " %d , %n", &numbers_of_residues_to_be_read[i], &N_of_characters_read);
                    delta = delta + N_of_characters_read;
                    
                    if (numbers_of_residues_to_be_read[i] == 0)
                        break;
                    total_number_of_residues_to_be_read++;
                }
            }
            else
            {
                printf("Something wrong with reading a configuration file\n");
            }
        }
        printf("Total number of residues to be read is %d\n", total_number_of_residues_to_be_read);
        printf("%d %d %d %d %d\n", numbers_of_residues_to_be_read[0], numbers_of_residues_to_be_read[1], numbers_of_residues_to_be_read[2], numbers_of_residues_to_be_read[3], numbers_of_residues_to_be_read[4]);
        if (fclose(fc) == EOF)
        {
            printf("Error closing the file %s\n", input);
            if (errno != 0)
                printf("Explanation: %s\n", strerror(errno));
            return 1;
        }
        return 0;
}

void draw_graph()
{
    int i = 0;
    int j = 0;
    int k = 0;
    int n = 0;
    int m = 0;
    //Width of legend field + some space between the legend and the graph
    int l_width = 0;
    //Width of one segment of histogram
    int s_width = 0;
    int border = 0;
    int gap = 2;
    int x_start_graph = 0;
    int x_start_labeling = 0;
    int x_start_histogram = 0;
    int y_start = 0;
    int y_start_labeling = 0;
    int x_start_legend = 0;
    int y_start_legend = 0;
    int label_x_number = 0;
    char label_x_text[5] = "";
    double scaling = 1.25;
    char PDB_filename[MAXSIZE] = "";


    dev2 = g2_open_X11(width_from_configuration_file, height_from_configuration_file);
    l_width = width_from_configuration_file/8;
    if (width_from_configuration_file > height_from_configuration_file)
    {
        scaling = (width_from_configuration_file/height_from_configuration_file) * 1.75;
    }
    border = 40 / (0.75 * scaling);
    s_width = ((width_from_configuration_file - (2 * border) - (100 * gap) - l_width)/scaling)/50;
    x_start_graph = border + (((s_width + 2 * gap)/2) * (scaling/1.25));
    x_start_histogram = border + gap;
    y_start = height_from_configuration_file - height_from_configuration_file/10 - border;
    y_start_labeling = y_start - ((s_width + (2 * gap)) * (scaling/1.25));
    x_start_labeling = border;
    x_start_legend = width_from_configuration_file - (0.75 * l_width);
    y_start_legend = y_start + (height_from_configuration_file/15);

    g2_set_line_width(dev2, 1);
    //Loop for drawing various graph elements: background, one-letter codes, graph of accessability, histogram etc.
    for (i = 0; i < (N_residues/50.0); i++)
    {
        for (j = n; ((j < (n + 50)) && (j < (N_residues - 1))); j++, m++)
        {
            //Drawing gray background rectangles
            g2_set_line_width(dev2, 1);
            g2_pen(dev2, g2_ink(dev2, 220/255.0, 220/255.0, 220/255.0));
            g2_filled_rectangle(dev2, x_start_labeling + (((s_width + (2 * gap)) * m) * (scaling/1.25)) - (gap/2), y_start - k, x_start_labeling + (((s_width + (2 * gap)) * (m + 1)) * (scaling/1.25)) + (gap/2), y_start - k + (height_from_configuration_file/15));
            //Drawing rectangles for letters
            g2_set_line_width(dev2, 1);
            g2_pen(dev2, g2_ink(dev2, residue_types[residues[j].residue_type].color_r, residue_types[residues[j].residue_type].color_g, residue_types[residues[j].residue_type].color_b));
            g2_filled_rectangle(dev2, x_start_labeling + (((s_width + (2 * gap)) * m) * (scaling/1.25)), y_start_labeling - k, x_start_labeling + (((s_width + (2 * gap)) * (m + 1)) * (scaling/1.25)), y_start - k);
            //Adding one-letter codes
            //Dealing with letter contrast. For dark colors of rectangles, letter color should be white, otherwise - black
            if (residues[j].residue_type == 2 || residues[j].residue_type == 6 || residues[j].residue_type == 10 || residues[j].residue_type == 11 || residues[j].residue_type == 20)
            {
                g2_pen(dev2, 0);
            }
            else
            {
                g2_pen(dev2, 1);
            }
            g2_set_font_size(dev2, ((s_width + (2 * gap))));
            g2_string(dev2, x_start_labeling + ((width_from_configuration_file - l_width)/500) + (((s_width + (2 * gap)) * m) * (scaling/1.25)), y_start_labeling + (s_width/4) - k, &residue_types[residues[j].residue_type].code1);
            //Drawing histogram of hydrophobicity
            g2_set_line_width(dev2, 1);
            g2_pen(dev2, g2_ink(dev2, 125/255.0, 127/255.0, 250/255.0));
            g2_filled_rectangle(dev2, x_start_histogram + (((s_width + (2 * gap)) * m) * (scaling/1.25)), y_start - k, x_start_histogram + (((s_width + (2 * gap)) * m) * (scaling/1.25)) + ((s_width + (gap)) * (scaling/1.25)), y_start - k + ((-(residue_types[residues[j].residue_type].hydrophobicity) + 5) * (height_from_configuration_file/150)));
            //Drawing x axis labels
            if (m == 0 || ((m + 1) % 10) == 0)
            {
                g2_pen(dev2, 1);
                label_x_number = j + 1;
                sprintf(label_x_text, "%d", label_x_number);
                g2_string(dev2, x_start_labeling + gap + ((width_from_configuration_file - l_width)/500) + (((s_width + (2 * gap)) * m) * (scaling/1.25)), y_start_labeling - ((height_from_configuration_file)/40) - k, label_x_text);
            }
            //Drawing x axis marking
            if (m == 0 || ((m + 1) % 5) == 0)
            {
                g2_set_line_width(dev2, 1);
                g2_pen(dev2, 1);
                g2_line(dev2, x_start_labeling + (((s_width + (2 * gap)) * (m + 0.5)) * (scaling/1.25)), y_start_labeling - k, x_start_labeling + (((s_width + (2 * gap)) * (m + 0.5)) * (scaling/1.25)), y_start_labeling - (((height_from_configuration_file) * (scaling/3.5))/100) - k);
            }
            //Drawing graph of accessability
            g2_pen(dev2, 19);
            g2_set_line_width(dev2, 3);
            //Since a histogram from the next cycle step will overlap some fragment of graph, it is necessary to either create another cycle for graph drawing or to implement some kind of delay in graph drawing. Probably, there are some other options, but I fount only two of them
            if (m == 0)
                continue;
            g2_line(dev2, x_start_graph + (((s_width + (2 * gap)) * (m - 1)) * (scaling/1.25)), y_start - k + (double)((double)height_from_configuration_file/12.5) - ((double)residues[j - 1].accessability * (double)((double)height_from_configuration_file/1000.0)), x_start_graph + (((s_width + (2 * gap)) * m) * (scaling/1.25)), y_start - k + (double)((double)height_from_configuration_file/12.5) - ((double)residues[j].accessability * (double)((double)height_from_configuration_file/1000.0)));
            g2_set_line_width(dev2, 1);
            
        }
        //Drawing black borders for grey background rectangles
        g2_pen(dev2, 1);
        g2_rectangle(dev2, x_start_labeling - (gap/2), y_start - k, x_start_labeling + (((s_width + (2 * gap)) * m) * (scaling/1.25)) + (gap/2), y_start - k + (height_from_configuration_file/15));
        
        //Drawing min and max labels
        g2_pen(dev2, 1);
        g2_set_font_size(dev2, ((s_width) + (2 * gap)));
        g2_string(dev2, x_start_graph - (2 * s_width) - (7 * gap), y_start - k - (s_width/5), "min");
        g2_string(dev2, x_start_graph - (2 * s_width) - (7 * gap), y_start - k - (s_width/5) + (height_from_configuration_file/15), "max");
        //Drawing hydrophobic and hydrophilic labels
        g2_string(dev2, x_start_graph + (((s_width + (2 * gap)) * m) * (scaling/1.25)), y_start - k - (s_width/5), "hydrophobic");
        g2_string(dev2, x_start_graph + (((s_width + (2 * gap)) * m) * (scaling/1.25)), y_start - k - (s_width/5) + (height_from_configuration_file/15), "hydrophilic");
        
        //Modifying variables, which control cycle and its parameters
        k = k + ((double)height_from_configuration_file/8.8) + s_width;
        n = n + 50;
        m = 0;
    }
    //Drawing legend
    i = 0;
    k = 0;
    for (i = 0; i <= 20; i++)
    {
        g2_pen(dev2, g2_ink(dev2, residue_types[i].color_r, residue_types[i].color_g, residue_types[i].color_b));
        g2_filled_rectangle(dev2, x_start_legend, y_start_legend - ((height_from_configuration_file/40) * i) - (height_from_configuration_file/100), x_start_legend + (width_from_configuration_file/50), y_start_legend - ((height_from_configuration_file/40) * (i + 1)));
        g2_set_font_size(dev2, ((s_width + (2 * gap))));
        g2_pen(dev2, 1);
        g2_string(dev2, x_start_legend + (width_from_configuration_file/35), y_start_legend - ((height_from_configuration_file/40) * (i + 0.5)) - (height_from_configuration_file/100), residue_types[i].code3);
    }
    
    //Drawing graph name and filename
    g2_set_font_size(dev2, (((2 * s_width) + gap)));
    g2_string(dev2, x_start_graph - (2 * gap), height_from_configuration_file - border, "Graph of hydrophobicity and accessability of residues");
    sprintf(PDB_filename, "PDB file: %s", input_to_reading);
    g2_string(dev2, x_start_graph - (2 * gap), border/3, PDB_filename);
}

void greeting()
{
    printf("\n");
    printf("--------------\n");
    printf("Final exercise\n");
    printf("--------------\n");
    printf("\n");
    printf("This program have the following features:\n1. ");
    printf("This program have the following features:\n1. ");
    printf("\n");
    printf("===========================================\n");
    printf("Honest warning, or better to say, confession\n");
    printf("===========================================\n");
    printf("Scaling has been implemented, but, to be honest, it is not perfect, so it would be better to test the program whithin the resolution range of [800; 1100] for both width and height (e.g. 800 x 1100, 1100 x 800, 800 x 800, 1100 x 1100, or any other value between 800 and 1100). The reason for such constraints is that currently the program in its scaling part partially rests on empirical values, and only partially on some variables, derived from width and height values, extracted from a configuration file. I have an idea of some clever scaling function or algorithm here, but it would require almost complete rewriting of significant part of graph drawing function.\n");
    printf("\n");
    printf("============\n");
    printf("Notification\n");
    printf("============\n");
    printf("If you are sure that you have already supplied enough arguments (name of the configuration file) via the command line, press Enter to proceed. Otherwise, press Enter, and then you will be required to enter name of the configuration file\n");
    getchar();
}
        
        
        

int main(int argc, char *argv[])
{
    //It was necessary to use argc == 2, not argc > 1, since I included some kind of protection in case too many arguments (as well as not enough arguments) were supplied via the command line
    if (argc == 2)
    {
        greeting();
       
        int function3 = configuration_file_reading(argv[1]);
        int function1 = file_opening_and_reading(input_to_reading);
        first_atom_scout();
        barycentre_calc();
        backbone_atoms();
        C_alpha_scout();
        C_alpha_residues_comparing();
        int function2 = file_opening_and_writing(DEFAULT_OUTPUT);
        draw_graph();
        getchar();
        g2_close(dev2);
        if (function1 == 1)
        {
            return 1;
        }
        if (function2 == 1)
        {
            return 1;
        }
        if (function3 == 1)
        {
            return 1;
        }
    }
    //Protection against less than two arguments supplied via the command line
    //In case there would be no arguments, program still can run, but it will ask the user to specify input and output filenames
    else if (argc < 2)
    {
        greeting();
        
        printf("Not enough arguments supplied. Please, provide the name of the configuration file\n");
        char input_filename[MAXSIZE] = "";
        printf("Please, enter the name of the configuration file: ");
        scanf("%s", input_filename);
        getchar();
        
        int function3 = configuration_file_reading(input_filename);
        int function1 = file_opening_and_reading(input_to_reading);
        first_atom_scout();
        barycentre_calc();
        backbone_atoms();
        C_alpha_scout();
        C_alpha_residues_comparing();
        int function2 = file_opening_and_writing(DEFAULT_OUTPUT);
        draw_graph();
        getchar();
        g2_close(dev2);
        if (function1 == 1)
        {
            return 1;
        }
        if (function2 == 1)
        {
            return 1;
        }
        if (function3 == 1)
        {
            return 1;
        }
    }
    //Protection against more than two arguments supplied via the command line
    //In case there would be no arguments, program still can run, but it will ask the user to specify input and output filenames
    else
    {
        greeting();
        
        printf("Too many arguments supplied. Please, provide the name of the configuration file only\n");
        char input_filename[MAXSIZE] = "";
        printf("Please, the name of the configuration file: ");
        scanf("%s", input_filename);
        getchar();
        
        int function3 = configuration_file_reading(input_filename);
        int function1 = file_opening_and_reading(input_to_reading);
        first_atom_scout();
        barycentre_calc();
        backbone_atoms();
        C_alpha_scout();
        C_alpha_residues_comparing();
        int function2 = file_opening_and_writing(DEFAULT_OUTPUT);
        draw_graph();
        getchar();
        g2_close(dev2);
        if (function1 == 1)
        {
            return 1;
        }
        if (function2 == 1)
        {
            return 1;
        }
        if (function3 == 1)
        {
            return 1;
        }
    }
    printf("\n");
    printf("-------------------------------\n");
    printf("Thank you for using the program\n");
    printf("-------------------------------\n");
    printf("\n");
    return 0;
}

    
    
    

/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0
*/

/*
  An FH5 to tecplot converter. This only works for specially written
  FH5 files. (This is designed primarily to work with TACS).
*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>

// // Include FH5 header files
// #include "FH5.h"

// Include FH5 header files
#include "TACSFH5Loader.h"
#include "TACSElementTypes.h"


// Include Tecplot header files
#include "TECIO.h"

enum FileType { FULL=0, GRID=1, SOLUTION=2 };

enum ZoneType { ORDERED=0, FELINESEG, FETRIANGLE, 
                FEQUADRILATERAL, FETETRAHEDRON, FEBRICK, 
                FEPOLYGON, FEPOLYHEDRA };

/*
  Initialize a data file.

  data_info == A description of the data
  var_names == Comma separated variable names
  file_name == Name of the file
  dir_name  == Name of the directory
  file_type == Type of data
*/
int create_tec_file( char *data_info, char *var_names,
                     char *file_name, char *dir_name,
                     enum FileType _file_type ){
  INTEGER4 file_type = _file_type;
  INTEGER4 debug = 0;
  INTEGER4 variables_are_double = 0;
  return TECINI112(data_info, var_names, file_name, dir_name,
                   &file_type, &debug, &variables_are_double);
}

/*
  A method to create a zone without the extra stuff
  
  zone_name    == The name of the zone to use
  zone_type    == One of LINESEG, TRIANGLE, QUAD, BRICK etc
  num_points   == The number of points
  num_elements == The number of elements
*/
int create_fe_tec_zone( char *zone_name, ZoneType _zone_type,
                        int _num_points, int _num_elements,
                        int use_strands=0,
                        double solution_time=0.0 ){
  if (_zone_type == ORDERED ||
      _zone_type == FEPOLYGON ||
      _zone_type == FEPOLYHEDRA ){
    fprintf(stderr, "Cannot create finite element zone with given \
zone type\n");
    return -1;
  }

  INTEGER4 zone_type = _zone_type;
  INTEGER4 num_points = _num_points;
  INTEGER4 num_elements = _num_elements;
  INTEGER4 num_faces = 0; // For all zones allowed here
  INTEGER4 icmax = 0, jcmax = 0, kcmax = 0; // Ignored
  INTEGER4 strand_id = 0;
  INTEGER4 parent_zone = 0;
  INTEGER4 is_block = 1; // Apparently this always needs to be 1
  // These are only for cell-based finite element data - we use node-based
  INTEGER4 num_face_connections = 0;
  INTEGER4 face_neighbour_mode = 0;
  INTEGER4 total_num_face_nodes = 0;
  INTEGER4 num_connected_boundary_faces = 0;
  INTEGER4 total_num_boundary_connections = 0;
  INTEGER4 *passive_var_list = NULL; // No passive variables
  INTEGER4 *value_location = NULL; // All values are nodal values
  INTEGER4 *share_var_from_zone = NULL;
  INTEGER4 share_con_from_zone = 0;

  // If we're using strands, set the strand ID
  if (use_strands){
    strand_id = 1;
  }

  return TECZNE112(zone_name, &zone_type, &num_points, &num_elements, 
                   &num_faces, &icmax, &jcmax, &kcmax,
                   &solution_time, &strand_id, &parent_zone, &is_block,
                   &num_face_connections, &face_neighbour_mode,
                   &total_num_face_nodes, &num_connected_boundary_faces,
                   &total_num_boundary_connections, 
                   passive_var_list, value_location,
                   share_var_from_zone, &share_con_from_zone);
}

/*
  Write data to a tecplot file

  len  == Length of the data
  data == The array of data
*/
int write_tec_double_data( int _len, double *data ){
  INTEGER4 len = _len;
  INTEGER4 is_double = 1;
  return TECDAT112(&len, data, &is_double);
}

/*
  Write float data to a tecplot file
*/
int write_tec_float_data( int _len, float *data ){
  INTEGER4 len = _len;
  INTEGER4 is_double = 0;
  return TECDAT112(&len, data, &is_double);
}

/*
  Write the connectivity data
*/
int write_con_data( int * con_data ){
  return TECNOD112(con_data);
}

/*
  End the file output
*/
int close_tec_file(){
  return TECEND112();
}

int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);

  // Check if we're going to use strands or not
  int use_strands = 0;
  int time_from_name = 0;
  
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "--use_strands") == 0){
      use_strands = 1;
    }
    if (strcmp(argv[k], "--time_from_name") == 0){
      time_from_name = 1;
    }
  }

  // Convert hdf5 file argv[1] to tecplot file argv[2]
  if (argc == 1){
    fprintf(stderr, "Error, no input files\n");
    return (1);
  }
  
  for ( int k = 1; k < argc; k++ ){
    char *infile = new char[ strlen(argv[k])+1 ];
    strcpy(infile, argv[k]);

    FILE *fp = fopen(infile, "r");
    if (fp){
      fclose(fp);
    }
    else {
      delete [] infile;
      continue;
    }

    // Set the output file
    char *outfile = new char[ strlen(infile)+5 ];
    int len = strlen(infile);
    int i = len-1;
    for ( ; i >= 0; i-- ){
      if (infile[i] == '.'){ break; }     
    }
    strcpy(outfile, infile);
    strcpy(&outfile[i], ".plt");
    
    printf("*NEW* Trying to convert FH5 file %s to tecplot file %s\n", 
           infile, outfile);

    char data_info[] = "Created by f5totec";
    char dir_name[] = "."; // We'll be working with the current directory
    int tec_init = 0;      // Tecplot initialized flag

    // // Open the FH5 file for reading
    // FH5File * file = new FH5File(MPI_COMM_SELF);
    // file->incref();
    // if (!file->openFile(infile)){
    //   fprintf(stderr, "Failed to open the file %s\n", infile);
    //   return (1);
    // }
    
        // Create the loader object
    TACSFH5Loader *loader = new TACSFH5Loader();
    loader->incref();


    int fail = loader->loadData(infile);
    if (fail){
      fprintf(stderr, "Failed to open the file %s\n", infile);
      return (1);
    }

      // Retrieve all the data from the file including the 
    // variables, connectivity and component numbers
    int num_elements;
    int *comp_nums, *ltypes, *ptr, *conn;
    loader->getConnectivity(&num_elements, &comp_nums, &ltypes, &ptr, &conn);

    const char *cname, *cvars;
    int num_nodes_continuous, num_vals_continuous;
    float *cdata;
    loader->getContinuousData(&cname, &cvars, &num_nodes_continuous, &num_vals_continuous, &cdata);


    const char *ename, *evars;
    int num_nodes_elements, num_vals_elements;
    float *edata;
    loader->getElementData(&ename, &evars, &num_nodes_elements, &num_vals_elements, &edata);



    int num_comp = loader->getNumComponents();
    int num_points = num_nodes_continuous;
    int conn_dim = 8;
    

    double *data = NULL;
    int *element_comp_num = comp_nums;
    
    
    // combine the continus variables and the element based variables
    int num_variables = num_vals_continuous + num_vals_elements ;
    float *float_data = new float[ num_variables*num_points ];
    memset(float_data, 0, num_variables*num_points*sizeof(float));
    
    
    for (int v = 0; v < num_vals_continuous; v++){
      for ( int k = 0; k < num_points; k++ ){
        float_data[v + k*num_variables] = cdata[v + k*num_vals_continuous];
      }
    }
    
    float *counts = new float[ num_points ];
    memset(counts, 0, num_points*sizeof(float));
    for ( int j = 0; j < ptr[num_elements]; j++ ){
      counts[conn[j]] += 1.0;
    }
    for ( int i = 0; i < num_points; i++ ){
      if (counts[i] != 0.0){
        counts[i] = 1.0/counts[i];
      }
    }
    
    //  For each component, average the nodal data
    for ( int v = num_vals_continuous; v < num_variables; v++ ){
      // Nodally average the data
      int j = v - num_vals_continuous;
      for ( int k = 0; k < ptr[num_elements]; k++ ){
        // data[conn[k]] += counts[conn[k]]*edata[edim2*k + j];
        // fprintf(stderr, "v:%d pt:%d , %d  idx:%d  edata: j:%d k:%d idx:%d", v, conn[k],conn[k]*num_variables, j, k, k*num_vals_elements + j );
        float_data[v + conn[k]*num_variables] += counts[conn[k]]*edata[k*num_vals_elements + j];
      }
    }
    

    double solution_time = 0.0;
    
    if (time_from_name){
      
      
      // determine the index of just the filename
      int j = strlen(infile);
      for ( ; j > 0; j-- ){
        if (infile[j] == '/'){ break; }     
      }
      char *filename = new char[i - j];
      
      // copy just the name of the file to a new string 
      strcpy(filename , &infile[j]);
      int k = 0;
      for (; j < i; j++){
        
        filename[k] = infile[j];
        k++;
      }
      
      float time; 
      printf(filename);
      printf("\n");
      //from stack overflow https://stackoverflow.com/a/4073314/14227912
      char *s = filename; 
  
      while (*s && !isdigit(*s)) s++; 
      
      if (*s)
      {
        sscanf(s, "%f", &time);  
        printf("time %f\n", time); 

      }    
        
      // cast the time taken from the file and set solution_time
      solution_time = (double)time;
      printf("solution_time %f\n", solution_time); 

    }
    
    //  Initialize the tecplot file with the variables
    char *vars = new char[ strlen(cvars)+ 2 + strlen(evars)+1 ];
    
    // this is the only way i could figure out how to concatenate these strings...
    strcpy(&vars[0], cvars);
    vars[strlen(cvars)] = ',';
    vars[strlen(cvars)+1] = ' ';
    strcpy(&vars[strlen(cvars)+ 2], evars);
    
    create_tec_file(data_info, vars,
                outfile, dir_name, FULL);
    tec_init = 1;
    delete [] vars;

  
    // Set the element type to use
    ZoneType zone_type;
    int convert[conn_dim];
    if (conn_dim == 2){      
      zone_type = FELINESEG; 
      convert[0] = 0;
      convert[1] = 1;
  }
    else if (conn_dim == 8){ 
      zone_type = FEBRICK; 
      convert[0] = 0;
      convert[1]  = 1;
      convert[2]  = 3;
      convert[3]  = 2;
      convert[4]  = 4;
      convert[5]  = 5;
      convert[6]  = 7;
      convert[7]  = 6;

    }
    else {
      zone_type = FEQUADRILATERAL; 
      convert[0] = 0;
      convert[1] = 1;
      convert[2] = 3;
      convert[3] =  2;

    }

    // int num_comp = file->getNumComponents();

    
    int *reduced_points = new int[ num_points ];
    int *reduced_conn = new int[ conn_dim*num_elements ];
    double *reduced_data = NULL;
    float *reduced_float_data = NULL;
    if (data){ 
      reduced_data = new double[ num_points ];
    }
    else if (float_data){
      reduced_float_data = new float[ num_points ];      
    }

    // fprintf(stderr, "num_elements:%d conn_dim:%d num_comp:%d num_points:%d\n",num_elements, conn_dim, num_comp, num_points );

    for ( int k = 0; k < num_comp; k++ ){
      // Count up the number of elements that use the connectivity
      char *comp_name = loader->getComponentName(k);
      //printf("Converting zone %d: %s at time %g\n", 
      // k, comp_name, solution_time);

      // fprintf(stderr, "k:%d comp_name:%s\n", k, comp_name);
    
      memset(reduced_points, 0, num_points*sizeof(int));
      memset(reduced_conn, 0, conn_dim*num_elements*sizeof(int));
      
      // fprintf(stderr, "reduced_conn = [\n");
      int npts = 1, nelems = 0;
      // Count up the number of points/elements in this sub-domain
      for ( int i = 0; i < num_elements; i++ ){
        if (element_comp_num[i] == k){
          // Add this element to the reduced connectivity
          // fprintf(stderr, "\n");
          for ( int j = 0; j < conn_dim; j++ ){
            int pt = conn[conn_dim*i + j];
            
          
            // If a reduced numbering has not been applied to this point,
            // create a new number for it
            // fprintf(stderr, "i:%d j:%d ", i,j);
            // fprintf(stderr, "pt:%d  ", pt);
            if (reduced_points[pt] == 0){
              reduced_points[pt] = npts;
              npts++;
            }
          
            // Set the reduced connectivity
            reduced_conn[conn_dim*nelems + convert[j]] = reduced_points[pt];
            // fprintf(stderr,"npts:%d    r_conn:%d\n",  npts, reduced_conn[conn_dim*nelems + j] );
          }
          nelems++;
        }
      }
      // fprintf(stderr, "]\n");
      // Since we started at npts = 1, we have one more point
      // than the actual number of points.
      npts--;
      

      // fprintf(stderr, "reduced_conn = [\n");
      if (nelems > 0 && npts > 0){
        // Create the zone with the solution time
        create_fe_tec_zone(comp_name, zone_type, npts, nelems,
                           use_strands, solution_time);

        // Retrieve the data
        for ( int j = 0; j < num_variables; j++ ){
          if (reduced_data){
            for ( int i = 0; i < num_points; i++ ){
              if (reduced_points[i] > 0){
                reduced_data[reduced_points[i]-1] = data[i*num_variables + j];
              }
            }
            write_tec_double_data(npts, reduced_data);
          }
          else if (reduced_float_data){
            for ( int i = 0; i < num_points; i++ ){
              if (reduced_points[i] > 0){
                reduced_float_data[reduced_points[i]-1] = 
                  float_data[i*num_variables + j];


              }
            }
            write_tec_float_data(npts, reduced_float_data);
          }
        }
        
        
        // Now, write the connectivity
        write_con_data(reduced_conn);
      }
    }

    if (tec_init){
      close_tec_file();
    }
    // file->close();
    // file->decref();
    loader->decref();

    // Clean up memory
    delete [] reduced_points;
    delete [] reduced_conn;
    delete [] reduced_data;

    delete [] data;

    // delete [] conn;
    // delete [] element_comp_num;

    delete [] infile;
    delete [] outfile;
  }
  MPI_Finalize();

  return (0);
}


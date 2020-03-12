/* one dimensional program for advection+dispersion, tested against ex11 in phreeqc
alpha = 1, upwind scheme for v > 0
alpha = 0, upwind scheme for v < 0
alpha = 0.5, central in space
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <IPhreeqc.h>
#include <mpi.h>
#include <math.h>
#define Nel 7  // number of transported element + H + O + chargebalance
#define Ncel 80 // solution number (40) 
#define MIN(a,b) ((a) < (b) ? (a) : (b))
int id;
int vt[Nel];
//double dValue[7][N];
double dv[Nel], cv[Nel];  //dv totoal mol, cv concentration
char dsv[Nel][100];
char sv[Nel][100];
char buffer[100];
//double por[7];
//double diff[7];
//double r_in[7];
//double r_out[7];
//double area_in[7];
//double area_out[7];
//double vol[7];
//double mix[7][2];
//bool radial;
//radial = 0;  // 0 for cartisian coordinate and 1 for radial coordinate 

void Extract(int cell)
{

                VAR v;
                int jj;
                VarInit(&v);

                        for(jj = 0; jj < Nel; jj++)

                                {

                                GetSelectedOutputValue(id,1,jj,&v);
                                vt[jj] = v.type;

                                switch (vt[jj])
                                                        {
                                                        case TT_DOUBLE:
                                                        dv[jj] = v.dVal;
                                                        sprintf(sv[jj],"%25.17e", v.dVal);
                                                        break;
                                                        case TT_STRING:
                                                        strcpy(sv[jj],v.sVal);
                                                        break;
                                                        }
                }

                VarClear(&v);


}


char *ConCat(const char *str1, const char *str2)
{
	strcpy(buffer, str1);
	return strcat(buffer, str2);
}
void EHandler(void)
{
	OutputErrorString(id);
	exit(EXIT_FAILURE);
}

int main()
{
	double radius_out = 130; //vitrified HLW outer_radius of buffer 
	double radius_in = 0.0;//vitrified HLW outer_radius of overpack
	int outlet_bc = 1, inlet_bc = 1; // 1 const, 2 impossilble, 3 flux
	int coord_sys = 1; // 1 = cartesian, 2 = radial
	double tot_time = 4 * 86400 * 365; // total simulation time, 4 years
	double por_v = 15.0 / (86400 * 365); // pore velocity
	//double lenCel = (radius_out - radius_in)/nCel;
	double vol_whole;//, vol_temp;
	double fbc = 1.;
	double deltaT = (radius_out - radius_in) / (Ncel * por_v); // curant number: deltaX = deltaT * por_v to minimize numerical dispersion
	double length = 1.;
	double alpha = 0;

	double disp_corr = 1;
	if (inlet_bc == 3)
		disp_corr += (1. / Ncel);
	if (outlet_bc == 3)
		disp_corr += (1. / Ncel);

	//double total_vol = 0; 
	int nMix = tot_time / deltaT; // equivalent to Nshifts in Phreeqc
	int mynode, size, up_node, down_node;
	//double starttime,endtime;
	char sdeltaT[10];
	sprintf(sdeltaT,"%8.2e",deltaT);
	int np;
	FILE *outputfile, *errfile;
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Status status;
	//starttime = MPI_Wtime();
	if(Ncel % size == 0)
		np = Ncel/size;
	else
		{
			if (mynode == 0)
				fprintf(stderr,"To evenly decompose the domain, the cell number should be divisible by the processor number!\n");	
			exit(1);
		}

	// including two ghost cells
	double por[np+2];  
	double diff[np+2], disp[np+2];
	double r_in_p[np+2], r_in[Ncel+2];
	double r_out_p[np+2], r_out[Ncel+2];
	double area_in_p[np+2];
	double area_out_p[np+2];
	double vol_p[np+2];
	double mix_d_up[np][Nel]; // mix due to diffusion
	double mix_d_down[np][Nel]; // mix due to advection
	double mix_a_up, mix_a_down;
	double dValue[2][np+2][Nel]; //first dimension for old and new solution
	char name[32];
	sprintf(name,"fig3_%d",mynode);

	//vol_whole = M_PI * (radius_out * radius_out - radius_in * radius_in) * length;

		for (int i = 0; i <= (Ncel+1); i++)   // assign outer and inner radii to each cell, calculate inner and outer surface areas 

			{
				if (i == 0)
				{
					r_in[i] = 0;
					r_out[i] = 0;
				}


				else if (i < (Ncel+1) )
				{
					r_in[i] = (i-1) * (radius_out - radius_in) / Ncel;
					r_out[i] = i * (radius_out - radius_in) / Ncel;
				}

				else
				{
					r_in[i] = r_out[i-1];
					r_out[i] = r_out[i-1];  //dummy zero width cell
				}
			}

		for (int i = 0; i <= (np+1); i++)   //initialize porosity and geometrical parameters. Here only two materials, TODO: automatically adapted

			 {
				r_out_p[i] = r_out[i+mynode*np]; 
				r_in_p[i] = r_in[i+mynode*np];
				area_out_p[i] =  M_PI * 0.04 * 0.04;  // diameter 0.08 m
				area_in_p[i] =  M_PI * 0.04 * 0.04;
				vol_p[i] =  M_PI * 0.04 * 0.04 * (r_out_p[i] - r_in_p[i]);
				vol_whole += vol_p[i];
		       }

	id = CreateIPhreeqc();
	SetOutputFileOn(id,1);
	SetLogFileOn(id,0);
	SetErrorFileOn(id,1);
	SetSelectedOutputFileOn(id,0);

	if (LoadDatabase(id, "phreeqc.dat") != 0) EHandler();
	if (RunFile(id, "fig2phreeqc2") != 0) EHandler();
//initialize the cell
	double mixf; // b_factor is boundary factor for flux boundary (3) compared to const boundary (1) 
	int dMix; 
	double diff_l; // hydrodynamic dispersion coefficient

	double bc_in_factor, bc_out_factor;
	

	diff_l = 0 * 2e-9 + 5 * por_v; 
	mixf = disp_corr * diff_l * deltaT / pow((radius_out - radius_in) / Ncel, 2);
	dMix = 1 + trunc(mixf * 3);
	mixf = mixf / dMix; 

	if (inlet_bc == 3)
		bc_in_factor = 0; // (por_v * (radius_out - radius_in) / Ncel) / (por_v * (radius_out - radius_in) / Ncel + diff_l * 2);
	else
		bc_in_factor = 1;
	if (outlet_bc == 3)
		bc_out_factor = 0; // (por_v * (radius_out - radius_in) / Ncel) / (por_v * (radius_out - radius_in) / Ncel + diff_l * 2);
	else
		bc_out_factor = 1;

	for (int j=0; j<= np+1; j++) 
	{
		char cj[10];
		sprintf(cj,"%d",j+mynode*np);
		AccumulateLine(id, ConCat("RUN_CELLS; -cells;",cj));
		AccumulateLine(id,";END");
		if (RunAccumulated(id) != 0) EHandler();
			Extract(j+mynode*np);

		for(int i=0; i < Nel; i++) 

			{
				dValue[0][j][i] = dv[i];  // initialize cell concentration

				if((mynode * np + i) == 0 || (mynode * np + i ) == (Ncel + 1))
					{
						mix_d_up[j][i] = mixf; 
						mix_d_down[j][i] = mixf; 
					}
				else
					{
						mix_d_up[j][i] = mixf;  
						mix_d_down[j][i] = mixf; 
					}
			}

		disp[j] = 5; 
	}

	up_node = mynode + 1;
	if (up_node >= size) up_node = MPI_PROC_NULL;

	down_node = mynode - 1;
	if (down_node < 0) down_node = MPI_PROC_NULL;


	if ((mynode % 2) == 0)
		{
			MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 0, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 0, MPI_COMM_WORLD, &status); 

			MPI_Sendrecv(r_in_p + np, 1, MPI_DOUBLE, up_node, 22, r_out_p + (np+1), 1, MPI_DOUBLE, up_node, 22, MPI_COMM_WORLD, &status);
		}
	else
		{
			MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 0, dValue[0][0], Nel, MPI_DOUBLE, down_node, 0, MPI_COMM_WORLD, &status);

			MPI_Sendrecv(r_out_p + 1, 1, MPI_DOUBLE, down_node, 22, r_in_p, 1, MPI_DOUBLE, down_node, 22, MPI_COMM_WORLD, &status);
		}

	if ((mynode % 2) == 1)
		{
			MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 1, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 1, MPI_COMM_WORLD, &status); 

			MPI_Sendrecv(r_in_p + np, 1, MPI_DOUBLE, up_node, 33, r_out_p + (np+1), 1, MPI_DOUBLE, up_node, 33, MPI_COMM_WORLD, &status);
		}
	else
		{
			MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 1, dValue[0][0], Nel, MPI_DOUBLE, down_node, 1, MPI_COMM_WORLD, &status);

			MPI_Sendrecv(r_out_p + 1, 1, MPI_DOUBLE, down_node, 33, r_in_p, 1, MPI_DOUBLE, down_node, 33, MPI_COMM_WORLD, &status);
		}

MPI_Barrier(MPI_COMM_WORLD); //probably not important here.

for (int m=1; m <= nMix; m++)  //advective time update
	{
		for (int j=1; j<= np; j++) //initialize cells

			{	
				char cj[10];
				sprintf(cj,"%d",j+mynode*np);
			// advection starts from here
				if  (mynode == size - 1 && j == np) // outlet boundary condition; mix[][1] is mix with higher number cell, mix[][0] mix with low number cell
				{
					mix_a_down = 0;
					mix_a_up = 2 * por_v * deltaT * (1 - alpha) / (r_out_p[j] - r_in_p[j-1]);
				}
				else if (mynode == 0 && j == 1) // 1 cell from inlet
				{
					mix_a_down = (-2) * por_v * deltaT * alpha / (r_out_p[j+1] - r_in_p[j]);

					mix_a_up =  1 * por_v * deltaT * (1 - alpha) / (r_out_p[j] - r_in_p[j]);// to make sure  por_v * deltaT = deltaX;

				}
				
				else
				{
					mix_a_up = 2 * por_v * deltaT * (1 - alpha) / (r_out_p[j] - r_in_p[j-1]);
					mix_a_down = (-2) * por_v * deltaT * alpha / (r_out_p[j+1] - r_in_p[j]);
				}
			//	printf("mix_a_up and mix_a_down are %f and %f\n",mix_a_up, mix_a_down);
				for (int k=0; k< Nel; k++) 
					{
					       dValue[1][j][k] = mix_a_down * dValue[0][j+1][k] + (1 - mix_a_down - mix_a_up) * dValue[0][j][k] + mix_a_up * dValue[0][j-1][k]; 

					       if(k != 2 && dValue[1][j][k] <= -0)
							dValue[1][j][k] = 0;
					       sprintf(sv[k], "%23.15e", dValue[1][j][k]); 
					}	
			// reaction after advection
				AccumulateLine(id, ConCat("SOLUTION_MODIFY ", cj));
				AccumulateLine(id, ConCat("-total_h ", sv[0]));
				AccumulateLine(id, ConCat("-total_o ", sv[1]));
				AccumulateLine(id, ConCat("-cb ", sv[2]));
				AccumulateLine(id, ConCat("-totals ", ""));
				AccumulateLine(id, ConCat(" Na ", sv[3]));
				AccumulateLine(id, ConCat(" Cl ", sv[4]));
				AccumulateLine(id, ConCat(" K ", sv[5]));
				AccumulateLine(id, ConCat(" N(5) ", sv[6]));
				AccumulateLine(id, ConCat("RUN_CELLS; -cells;", cj));
				//AccumulateLine(id, "-time_step  1800");
				AccumulateLine(id, ";END");
				if (RunAccumulated(id) != 0) EHandler();
				Extract(j+mynode*np);
			// retardation factor calculated after each advection timestep, and before dispersive timestep 
				double Cdif[Nel], Rf[Nel];
				for (int k = 0; k < Nel; k++)
				{
			/*		if ( k < 3)
						Rf[k] = 1;  // no retardation for H, O, charge
					else
					{
						Cdif[k] = dv[k] - dValue[0][j][k]; 
						if (abs(Cdif[k] > 1e-9)) 
						//	Rf[k] = (dValue[1][j][k] - dValue[0][j][k]) / Cdif[k];
							Rf[k] = 1;
						if (Rf[k] < 1) 
							 Rf[k] = 1;
						// recalculate mixing factor considering the numerical dispersion caused by retardation.
					}
			*/
					mix_d_up[j][k] = mixf; // - (1 - por_v * deltaT / (Rf[k] * (radius_out - radius_in) / Ncel)) / (2 * dMix);
					if (mix_d_up[j][k] < 0) 
						mix_d_up[j][k] = 0;
					mix_d_down[j][k] = mixf; //- (1 - por_v * deltaT / (Rf[k] * (radius_out - radius_in) / Ncel)) / (2 * dMix);
					if (mix_d_down[j][k] < 0)
						mix_d_down[j][k] = 0;
					dValue[1][j][k] = dv[k]; // new concentration for cell j, not necessary here unless for output
				//	printf("for %ith cell, element %i, mixing factor is %f\n", j, k, mix_d_down[j][k]);
				}

	/*			if ((mynode*np + j) == Ncel) // for output results after advection
					{
						outputfile = fopen(name,"a");
						fprintf(outputfile,"%i %i %e %e %e %e %e %e %e %e %e\n",m,mynode*np+j,dv[0],dv[1],dv[2],dv[3],dv[4],dv[5],dv[6],dv[7],(m+0.5)/40); 
						fclose(outputfile);
					}   */

			}   // finish loop cells for one advetive timestep

			for (int j = 1 ; j <= np; j++)
				for (int k = 0; k < Nel; k++)
					dValue[0][j][k] = dValue[1][j][k];

	MPI_Barrier(MPI_COMM_WORLD); //probably not important here.

        // message passing is necessary here because dValue[][][] from neigbouring nodes are needed
	if ((mynode % 2) == 0)
		{
			MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 0, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 0, MPI_COMM_WORLD, &status); 

	//		MPI_Sendrecv(r_in_p + np, 1, MPI_DOUBLE, up_node, 22, r_out_p + (np+1), 1, MPI_DOUBLE, up_node, 22, MPI_COMM_WORLD, &status);
		}
	else
		{
			MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 0, dValue[0][0], Nel, MPI_DOUBLE, down_node, 0, MPI_COMM_WORLD, &status);

	//		MPI_Sendrecv(r_out_p + 1, 1, MPI_DOUBLE, down_node, 22, r_in_p, 1, MPI_DOUBLE, down_node, 22, MPI_COMM_WORLD, &status);
		}

	if ((mynode % 2) == 1)
		{
			MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 1, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 1, MPI_COMM_WORLD, &status); 

	//		MPI_Sendrecv(r_in_p + np, 1, MPI_DOUBLE, up_node, 33, r_out_p + (np+1), 1, MPI_DOUBLE, up_node, 33, MPI_COMM_WORLD, &status);
		}
	else
		{
			MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 1, dValue[0][0], Nel, MPI_DOUBLE, down_node, 1, MPI_COMM_WORLD, &status);

	//		MPI_Sendrecv(r_out_p + 1, 1, MPI_DOUBLE, down_node, 33, r_in_p, 1, MPI_DOUBLE, down_node, 33, MPI_COMM_WORLD, &status);
		}
			 
					// dispersive timestep, dMix is the number of dispersive mixing in each advection 

					for (int dm = 0; dm < dMix; dm++)
					    	{  
							for (int j=1; j<= np; j++) //initialize cells

								{	
									char kj[10];
									sprintf(kj,"%d",j+mynode*np);
									for (int k = 0; k < Nel; k++) 
									{
										if (mynode == 0 && j == 1)
											dValue[1][1][k] = mix_d_down[1][k] * dValue[0][2][k] + bc_in_factor *  mix_d_up[1][k] * dValue[0][0][k] + (1 - mix_d_down[1][k] - bc_in_factor * mix_d_up[1][k]) * dValue[0][1][k];
										else if (mynode == (size -1) && j == np)
											dValue[1][np][k] = mix_d_up[np][k] * dValue[0][np-1][k] + bc_out_factor * mix_d_down[np][k] * dValue[0][np+1][k] + (1 - mix_d_up[np][k] - bc_out_factor * mix_d_down[np][k]) * dValue[0][np][k];

										else
											dValue[1][j][k] = mix_d_down[j][k] * dValue[0][j+1][k]  + (1 - mix_d_up[j][k] - mix_d_down[j][k]) * dValue[0][j][k] + mix_d_up[j][k] * dValue[0][j-1][k]; 

										if ( k != 2 && dValue[1][j][k] <= -0)
											dValue[1][j][k] = 0;

										sprintf(dsv[k], "%23.15e", dValue[1][j][k]); 
									}

								// Reaction after dispersive transport

									AccumulateLine(id, ConCat("SOLUTION_MODIFY ", kj));
									AccumulateLine(id, ConCat("-total_h ", dsv[0]));
									AccumulateLine(id, ConCat("-total_o ", dsv[1]));
									AccumulateLine(id, ConCat("-cb ", dsv[2]));
									AccumulateLine(id, ConCat("-totals ", ""));
									AccumulateLine(id, ConCat(" Na ", dsv[3]));
									AccumulateLine(id, ConCat(" Cl ", dsv[4]));
									AccumulateLine(id, ConCat(" K ", dsv[5]));
									AccumulateLine(id, ConCat(" N(5) ", dsv[6]));
									AccumulateLine(id, ConCat("RUN_CELLS; -cells;", kj));
									//AccumulateLine(id, "-time_step  1800");
									AccumulateLine(id, ";END");
									if (RunAccumulated(id) != 0) EHandler();
										Extract(j+mynode*np);

									for (int k = 0; k < Nel; k++)
										{
											dValue[1][j][k] = dv[k];
										//	printf("cell %i,element %i,  value is %f\n", dj, k,dValue[0][dj][k]);
										}
									if (m == nMix && dm == dMix - 1)
										{
											outputfile = fopen(name,"a");
											fprintf(outputfile,"%i %i %f %e %e %e %e %e\n",m,mynode*np+j,(r_in_p[j] + r_out_p[j]) * 0.5,dv[0],dv[1],dv[2],dv[3],dv[4]); 
											fclose(outputfile);
										}
								} // loop over cell for dispersion + reaction

								// update new concentrations for all the cells
 
								for (int j = 1 ; j <= np; j++)
									for (int k = 0; k < Nel; k++)
										dValue[0][j][k] = dValue[1][j][k];

							MPI_Barrier(MPI_COMM_WORLD); //probably not important here.
							// mpi should be added here also.
							if ((mynode % 2) == 0)
								{
									MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 0, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 0, MPI_COMM_WORLD, &status); 

							//		MPI_Sendrecv(r_in_p + np, 1, MPI_DOUBLE, up_node, 22, r_out_p + (np+1), 1, MPI_DOUBLE, up_node, 22, MPI_COMM_WORLD, &status);
								}
							else
								{
									MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 0, dValue[0][0], Nel, MPI_DOUBLE, down_node, 0, MPI_COMM_WORLD, &status);

							//		MPI_Sendrecv(r_out_p + 1, 1, MPI_DOUBLE, down_node, 22, r_in_p, 1, MPI_DOUBLE, down_node, 22, MPI_COMM_WORLD, &status);
								}

							if ((mynode % 2) == 1)
								{
									MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 1, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 1, MPI_COMM_WORLD, &status); 

							//		MPI_Sendrecv(r_in_p + np, 1, MPI_DOUBLE, up_node, 33, r_out_p + (np+1), 1, MPI_DOUBLE, up_node, 33, MPI_COMM_WORLD, &status);
								}
							else
								{
									MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 1, dValue[0][0], Nel, MPI_DOUBLE, down_node, 1, MPI_COMM_WORLD, &status);

							//		MPI_Sendrecv(r_out_p + 1, 1, MPI_DOUBLE, down_node, 33, r_in_p, 1, MPI_DOUBLE, down_node, 33, MPI_COMM_WORLD, &status);
								}



						} // end of dMix times of dispersion + reactin 
		  // finish one step of advection + dMix dispersion

	MPI_Barrier(MPI_COMM_WORLD); //probably not important here.

	if ((mynode % 2) == 0)
	{
		MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 0, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 0, MPI_COMM_WORLD, &status); 

		MPI_Sendrecv(por + np, 1, MPI_DOUBLE, up_node, 222, por + (np+1), 1, MPI_DOUBLE, up_node, 222, MPI_COMM_WORLD, &status); //not important w no-update

		MPI_Sendrecv(diff + np, 1, MPI_DOUBLE, up_node, 555, diff + (np+1), 1, MPI_DOUBLE, up_node, 555, MPI_COMM_WORLD, &status);//not important wt update 
	}
	else
	{
		MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 0, dValue[0][0], Nel, MPI_DOUBLE, down_node, 0, MPI_COMM_WORLD, &status);

		MPI_Sendrecv(por + 1, 1, MPI_DOUBLE, down_node, 222, por + 0, 1, MPI_DOUBLE, down_node, 222, MPI_COMM_WORLD, &status); //not important without update


		MPI_Sendrecv(diff + 1, 1, MPI_DOUBLE, down_node, 555, diff + 0, 1, MPI_DOUBLE, down_node, 555, MPI_COMM_WORLD, &status); //not important without update
	}

	if ((mynode % 2) == 1)
	{
		MPI_Sendrecv(dValue[0][np], Nel, MPI_DOUBLE, up_node, 1, dValue[0][np+1], Nel, MPI_DOUBLE, up_node, 1, MPI_COMM_WORLD, &status); 

		MPI_Sendrecv(por + np, 1, MPI_DOUBLE, up_node, 333, por + (np + 1), 1, MPI_DOUBLE, up_node, 333, MPI_COMM_WORLD, &status);//not important wno upd 

		MPI_Sendrecv(diff + np, 1, MPI_DOUBLE, up_node, 444, diff + (np + 1), 1, MPI_DOUBLE, up_node, 444, MPI_COMM_WORLD, &status);//not important wno upd 

	}
	else
	{
		MPI_Sendrecv(dValue[0][1], Nel, MPI_DOUBLE, down_node, 1, dValue[0][0], Nel, MPI_DOUBLE, down_node, 1, MPI_COMM_WORLD, &status);

		MPI_Sendrecv(por + 1, 1, MPI_DOUBLE, down_node, 333, por + 0, 1, MPI_DOUBLE, down_node, 333, MPI_COMM_WORLD, &status);//not important wno upd

		MPI_Sendrecv(diff + 1, 1, MPI_DOUBLE, down_node, 444, diff + 0, 1, MPI_DOUBLE, down_node, 444, MPI_COMM_WORLD, &status);//not important wno upd
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
}   // end of looping over time
MPI_Finalize();
if(DestroyIPhreeqc(id) != IPQ_OK) EHandler();
exit(EXIT_SUCCESS);
return 0;
} // end of the main program

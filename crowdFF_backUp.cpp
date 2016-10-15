// Conway's game of Life Cellular Automaton

extern "C"{
#include <OpenCAL/cal2D.h>
#include <OpenCAL/cal2DIO.h>
#include <OpenCAL/cal2DRun.h>
#include <OpenCAL/cal2DUnsafe.h>
#include <OpenCAL-GL/calgl2D.h>
#include <OpenCAL-GL/calgl2DWindow.h>
}
#include <stdlib.h>
#include <vector>
#include <iostream>

#define CROWD_DOOR -1
#define MAX_STATICFF_VAL 50

// declare CA, substate and simulation objects
struct CALModel2D* crowdFF_model;
struct Substates {
	struct CALSubstate2Di* pedestrians; 
	struct calSaveSubstate2Db * conflicts; 
	struct CALSubstate2Dr* staticFF;
	 
	struct CALSubstate2Di* trajectories;  
	struct CALSubstate2Di* hotSpot;  
} Q;

struct CALRun2D* crowdFF_simulation;

void crowdFFInit(struct CALModel2D* crowdFF_model){
	
	// set the whole substate to 0
	calLoadSubstate2Di(crowdFF_model, Q.pedestrians, "./geometry.txt");
	calInitSubstate2Dr(crowdFF_model, Q.staticFF, MAX_STATICFF_VAL);
	calInitSubstate2Di(crowdFF_model, Q.trajectories, 0);
	calInitSubstate2Di(crowdFF_model, Q.hotSpot, 0);
	
	for (int i=1; i < crowdFF_model->rows-1; i++){
		for (int j=1; j < crowdFF_model->columns-1; j++){
			if (calGet2Di(crowdFF_model, Q.pedestrians, i, j) == -1){
				calSetCurrent2Dr(crowdFF_model, Q.staticFF, i,j, 0);
			}
		}
	}

	// brute force version of Static FF calculation
	bool calculation_pending = true;
	while (calculation_pending){
		calculation_pending = false;
		for (int i=1; i < crowdFF_model->rows-1; i++){
			for (int j=1; j < crowdFF_model->columns-1; j++){
				if (calGet2Di(crowdFF_model, Q.pedestrians, i, j) == 0){
					CALreal min = MAX_STATICFF_VAL;
					for (int n=1; n<crowdFF_model->sizeof_X; n++){
						CALreal nei_value = calGetX2Dr(crowdFF_model, Q.staticFF, i, j, n);
						if (min > nei_value ){
							min = nei_value;
						}
					}
					min +=1;
					if (calGet2Dr(crowdFF_model, Q.staticFF, i,j) -0.5> min){
						calSetCurrent2Dr(crowdFF_model, Q.staticFF, i,j, min);
						calculation_pending = true;
					}
				}
			}
		}
	}

	// static FF calculation 
	/*std::vector<CALCell2D> toCheck;
	for (int i=1; i < crowdFF_model->rows-1; i++){
		for (int j=1; j < crowdFF_model->columns-1; j++){
			if (calGet2Di(crowdFF_model, Q.pedestrians, i,j) == CROWD_DOOR){
				calSet2Dr(crowdFF_model, Q.staticFF, i, j, 0);
				for (int n=1; n < crowdFF_model->sizeof_X; n++){
					CALCell2D nei;
					
					nei.i = crowdFF_model->X[n].i;
					nei.j = crowdFF_model->X[n].j;
				//	toCheck.push_back()
				}
			}
		}
	}
	*/

	calSet2Di(crowdFF_model, Q.pedestrians, 6, 5, 2);
	calSet2Di(crowdFF_model, Q.pedestrians, 6, 4, 2);
//	calSet2Di(crowdFF_model, Q.pedestrians, 10,12, 2);
//	calSet2Di(crowdFF_model, Q.pedestrians, 6,12, 2);

	calUpdate2D(crowdFF_model);

	


}

// The pedestrians's transition function
void pedestrianTransitionFunction(struct CALModel2D* crowdFF_model, int i, int j)
{
	if(calGet2Di(crowdFF_model, Q.pedestrians, i, j) == 2){
		//extra slow down
		Sleep(100);

		std::cout << i  << " " << j << "\n";
		// move to function
		CALreal min = calGet2Dr(crowdFF_model, Q.staticFF, i, j);
		int desiredCell = 0;
		for (int n=1; n<crowdFF_model->sizeof_X; n++){
			CALreal nei_value = calGetX2Dr(crowdFF_model, Q.staticFF, i, j, n);

			if (min > nei_value  && calGetX2Di(crowdFF_model, Q.pedestrians, i, j, n)<1){
				min = nei_value;
				desiredCell = n;
			}
		}
		std::cout <<" miniumu" <<  min << "\n";
		if(desiredCell !=0){
			calSet2Di(crowdFF_model, Q.pedestrians, i,j,0);
			calSet2Di(crowdFF_model, Q.trajectories, i,j,1);
			calSet2Di(crowdFF_model, Q.hotSpot, i,j,calGet2Di(crowdFF_model, Q.hotSpot, i,j)+1);
			int desiredCellType = calGetX2Di(crowdFF_model, Q.pedestrians, i,j, desiredCell);
			if (desiredCellType == 0){			
				calSetX2Di(crowdFF_model, Q.pedestrians, i,j, desiredCell, 2);
			}
			if (desiredCellType == -1){			
				std::cout << "I want to go outside\n";
			}
			//calSet2Di(crowdFF_model, Q.pedestrians, crowdFF_model->X[desiredCell].i, crowdFF_model->X[desiredCell].j, 2);
		}
	}
}

void exitFunction()
{
	// save the Q substate to file
	calSaveSubstate2Di(crowdFF_model, Q.pedestrians, "./simulationState.txt");
	calSaveSubstate2Di(crowdFF_model, Q.trajectories, "./trajectories.txt");
	calSaveSubstate2Di(crowdFF_model, Q.hotSpot, "./hotSpot.txt");
	calSaveSubstate2Dr(crowdFF_model, Q.staticFF, "./staticFF.txt");
	// finalize simulation and CA objects
	calRunFinalize2D(crowdFF_simulation);
	calFinalize2D(crowdFF_model);
}

int main(int argc, char** argv)
{
	struct CALGLDrawModel2D* drawModel = NULL;
	struct CALGLDrawModel2D* drawModelStaticFF = NULL;

	struct CALGLDrawModel2D* drawModelTraj = NULL;
	struct CALGLDrawModel2D* drawModelHotSpot = NULL;
	atexit(exitFunction);


	// define of the crowdFF_model CA and crowdFF_simulation simulation objects
	crowdFF_model = calCADef2D(16, 16, CAL_MOORE_NEIGHBORHOOD_2D, CAL_SPACE_FLAT, CAL_NO_OPT);
	crowdFF_simulation = calRunDef2D(crowdFF_model, 1, 20, CAL_UPDATE_IMPLICIT);

	// add the Q substate to the crowdFF_model CA
	Q.pedestrians = calAddSubstate2Di(crowdFF_model);
	Q.staticFF = calAddSingleLayerSubstate2Dr(crowdFF_model);
	Q.trajectories = calAddSubstate2Di(crowdFF_model);
	Q.hotSpot = calAddSubstate2Di(crowdFF_model);

	// add transition function's elementary process
	calAddElementaryProcess2D(crowdFF_model, pedestrianTransitionFunction);

	

	// simulation run
	calRunAddInitFunc2D(crowdFF_simulation, crowdFFInit);
	calRunInitSimulation2D(crowdFF_simulation);


	//calRun2D(crowdFF_simulation);
	//system("PAUSE");
	// just a visualization !

	// Initialize the viewer
	calglInitViewer("crowdFF viewer", 1.0f, 400, 400, 40, 40, CAL_TRUE, 1);
	// model1 definition
	drawModel = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "crowdFF", crowdFF_model, crowdFF_simulation);
	// Add nodes'
	calglAdd2Di (drawModel, NULL, &Q.pedestrians, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_STATIC);
	calglAdd2Di (drawModel, Q.pedestrians, &Q.pedestrians, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_GREEN_SCALE, CALGL_DATA_TYPE_STATIC);

	drawModelStaticFF = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "staticFF", crowdFF_model, crowdFF_simulation);
	calglAdd2Dr (drawModelStaticFF, NULL, &Q.staticFF, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_STATIC);
	calglAdd2Dr (drawModelStaticFF, Q.staticFF, &Q.staticFF, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_GREEN_SCALE, CALGL_DATA_TYPE_STATIC);
	//calglInfoBar2Dr(drawModelStaticFF , Q.staticFF, "staticFF", CALGL_TYPE_INFO_USE_GREEN_SCALE , 20, 120, 300, 40);

	drawModelTraj = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "trajectories", crowdFF_model, crowdFF_simulation);
	calglAdd2Di (drawModelTraj, NULL, &Q.trajectories, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_STATIC);
	calglAdd2Di (drawModelTraj, Q.trajectories, &Q.trajectories, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_GREEN_SCALE, CALGL_DATA_TYPE_STATIC);
	
	drawModelHotSpot = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "hot spots", crowdFF_model, crowdFF_simulation);
	calglAdd2Di (drawModelHotSpot, NULL, &Q.hotSpot, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_STATIC);
	calglAdd2Di (drawModelHotSpot, Q.hotSpot, &Q.hotSpot, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_GREEN_SCALE, CALGL_DATA_TYPE_STATIC);
	


	calglSetDisplayStep(1);

	calglMainLoop2D(argc, argv);

	

	return 0;
}

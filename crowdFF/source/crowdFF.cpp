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
#define MAX_STATICFF_VAL 200

#define STATIC_PARAM 0.9
#define DYNAMIC_PARAM 0.9

#define DFF_NUM_PHEROMONES 2
#define DFF_DIFUSION_PROB 70
#define DFF_VAPORATION_PROB 10

// declare CA, substate and simulation objects
struct CALModel2D* crowdFF_model;
struct Substates {
	struct CALSubstate2Di* pedestrians; 
	struct CALSubstate2Db* conflicts; 

	struct CALSubstate2Dr* floorField;
	struct CALSubstate2Dr* staticFF;
	struct CALSubstate2Di* dynamicFF;
	 
	struct CALSubstate2Di* trajectories;  
	struct CALSubstate2Di* hotSpot;  
} Q;

struct CALRun2D* crowdFF_simulation;

void crowdFFInit(struct CALModel2D* crowdFF_model){
	
	calLoadSubstate2Di(crowdFF_model, Q.pedestrians, "./geometry3.txt");	
	calInitSubstate2Db(crowdFF_model, Q.conflicts, 0);
	calInitSubstate2Dr(crowdFF_model, Q.staticFF, MAX_STATICFF_VAL);
	calInitSubstate2Dr(crowdFF_model, Q.floorField, 0);
	calInitSubstate2Di(crowdFF_model, Q.dynamicFF, 0);
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

	/*calSet2Di(crowdFF_model, Q.pedestrians, 8, 5, 2);
	calSet2Di(crowdFF_model, Q.pedestrians, 8, 6, 2);
	calSet2Di(crowdFF_model, Q.pedestrians, 10,12, 2);
	calSet2Di(crowdFF_model, Q.pedestrians, 6,12, 2);
	calSet2Di(crowdFF_model, Q.pedestrians, 4,13, 2);
	*/
	calUpdate2D(crowdFF_model);
} 

void clearConflicts(struct CALModel2D* crowdFF_model){
	calInitSubstate2Db(crowdFF_model, Q.conflicts, 0);
}

unsigned char transformNeighbourToMask(int nei){
	switch(nei){
		case 1:
		return 1;

		case 2:
		return 2;

		case 3: 
		return 4;

		case 4:
		return 8;

		case 5:
		return 16;

		case 6: 
		return 32;

		case 7:
		return 64;

		case 8: 
		return 128;
	}
}

int resolveConflictFromTheMask(unsigned char mask){
	std::vector<int> pedestrians;
	if (mask & 1){
		pedestrians.push_back(1);
	}
	if (mask & 2){
		pedestrians.push_back(2);
	}
	if (mask & 4){
		pedestrians.push_back(3);
	}
	if (mask & 8){
		pedestrians.push_back(4);
	}
	if (mask & 16){
		pedestrians.push_back(5);
	}
	if (mask & 32){
		pedestrians.push_back(6);
	}
	if (mask & 64){
		pedestrians.push_back(7);
	}
	if (mask & 128){
		pedestrians.push_back(8);
	}

	if (pedestrians.size() ==1){
		return pedestrians.at(0);
	}
	else return pedestrians.at(rand()%pedestrians.size());	
}

int opositeNeighbour(int nei){
	switch(nei){
		case 1:
		return 4;

		case 2:
		return 3;

		case 3: 
		return 2;

		case 4:
		return 1;

		case 5:
		return 7;

		case 6: 
		return 8;

		case 7:
		return 5;

		case 8: 
		return 6;
	}
}

void totalFloorField(struct CALModel2D* crowdFF_model, int i, int j){
	float staticValue = calGet2Dr(crowdFF_model, Q.staticFF, i, j);
	int dynamicValue = calGet2Di(crowdFF_model, Q.dynamicFF, i, j);
	float floorFieldValue = staticValue*STATIC_PARAM - dynamicValue * DYNAMIC_PARAM;
	calSet2Dr(crowdFF_model, Q.floorField, i, j, floorFieldValue);
}

// The pedestrians's transition function 1
void pedestrianTransitionFunction(struct CALModel2D* crowdFF_model, int i, int j)
{
	
	if(calGet2Di(crowdFF_model, Q.pedestrians, i, j) == 2){
	//	CALreal min = calGet2Dr(crowdFF_model, Q.staticFF, i, j);
		CALreal min = calGet2Dr(crowdFF_model, Q.floorField, i, j);
		int desiredCell = 0;
		for (int n=1; n<crowdFF_model->sizeof_X; n++){
			CALreal nei_value = calGetX2Dr(crowdFF_model, Q.floorField, i, j, n);
			if (min > nei_value  && calGetX2Di(crowdFF_model, Q.pedestrians, i, j, n)<1){
				min = nei_value;
				desiredCell = n;
			}
		}
		if(desiredCell !=0){
			unsigned char mask = transformNeighbourToMask (desiredCell);
			unsigned char currMask = calGetX2Db(crowdFF_model, Q.conflicts, i,j, desiredCell);
			unsigned char newMask = currMask | mask;
			calSetX2Db(crowdFF_model, Q.conflicts, i,j, desiredCell, newMask);
		}
	}
}

// conflict resolution transition function
void conflictResolutionFunction(struct CALModel2D* crowdFF_model, int i, int j)
{
	unsigned char mask = calGet2Db(crowdFF_model, Q.conflicts, i, j);
	if(mask != 0){
		//extra slow down
	//	Sleep(1);

		
		int movingPedestrian = opositeNeighbour (resolveConflictFromTheMask (mask));
		
		calSetX2Di(crowdFF_model, Q.pedestrians, i,j,movingPedestrian, 0);
		calSetX2Di(crowdFF_model, Q.trajectories, i,j,movingPedestrian, 1);
		calSetX2Di(crowdFF_model, Q.hotSpot, i,j, movingPedestrian, calGetX2Di(crowdFF_model, Q.hotSpot, i,j, movingPedestrian)+1);
		

		int type = calGet2Di( crowdFF_model, Q.pedestrians, i,j);	
		if (type == 0){
			calSet2Di(crowdFF_model, Q.pedestrians, i,j,2);
			calSetX2Di(crowdFF_model, Q.dynamicFF, i,j, movingPedestrian, calGetX2Di(crowdFF_model, Q.dynamicFF, i,j, movingPedestrian)+DFF_NUM_PHEROMONES);
		}
		if (type == -1){
			std::cout <<"exit \n";
		}
	}
}

void updateDynamicFloorField(struct CALModel2D* crowdFF_model, int i, int j){
	int dynamicFFValue = calGet2Di(crowdFF_model, Q.dynamicFF, i, j);	
//	for (int x=0; x< dynamicFFValue; x++){
	if (dynamicFFValue > 0){
		int prob = rand()%100 +1;
		if(prob < DFF_VAPORATION_PROB){
			calSet2Di(crowdFF_model, Q.dynamicFF, i,j, calGetNext2Di(crowdFF_model, Q.dynamicFF, i, j) -1);
		}		
		else if (prob < DFF_VAPORATION_PROB + DFF_DIFUSION_PROB ){
			calSet2Di(crowdFF_model, Q.dynamicFF, i,j, calGetNext2Di(crowdFF_model, Q.dynamicFF, i, j)  -1);
			int neiNum = rand()%8+1;
			if (calGetX2Di(crowdFF_model, Q.pedestrians,i,j, neiNum) ==0 ){
				int pheromoneNum = calGetX2Di(crowdFF_model, Q.dynamicFF,i,j, neiNum) +1;
				
				calSetX2Di(crowdFF_model, Q.dynamicFF,i,j, neiNum,  pheromoneNum);
			}
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
	calSaveSubstate2Di(crowdFF_model, Q.dynamicFF, "./dynamicFF.txt");
	// finalize simulation and CA objects
	calRunFinalize2D(crowdFF_simulation);
	calFinalize2D(crowdFF_model);
}

int main(int argc, char** argv)
{
	srand (time(NULL));

	struct CALGLDrawModel2D* drawModel = NULL;

	struct CALGLDrawModel2D* drawModelFloorField = NULL;
	struct CALGLDrawModel2D* drawModelStaticFF = NULL;
	struct CALGLDrawModel2D* drawModelDynamicFF = NULL;
	
	struct CALGLDrawModel2D* drawModelTraj = NULL;
	struct CALGLDrawModel2D* drawModelHotSpot = NULL;

	struct CALGLDrawModel2D* drawModelConflicts = NULL;
	atexit(exitFunction);


	// define of the crowdFF_model CA and crowdFF_simulation simulation objects
	crowdFF_model = calCADef2D(50, 50, CAL_MOORE_NEIGHBORHOOD_2D, CAL_SPACE_FLAT, CAL_NO_OPT);
	crowdFF_simulation = calRunDef2D(crowdFF_model, 1, 500, CAL_UPDATE_IMPLICIT);

	// add the Q substate to the crowdFF_model CA
	Q.pedestrians = calAddSubstate2Di(crowdFF_model);
	Q.conflicts = calAddSubstate2Db(crowdFF_model);

	Q.floorField = calAddSubstate2Dr(crowdFF_model);
	Q.staticFF = calAddSingleLayerSubstate2Dr(crowdFF_model);
	Q.dynamicFF = calAddSubstate2Di(crowdFF_model);

	Q.trajectories = calAddSubstate2Di(crowdFF_model);
	Q.hotSpot = calAddSubstate2Di(crowdFF_model);

	// add transition function's elementary process
	calAddElementaryProcess2D(crowdFF_model, totalFloorField);
	calAddElementaryProcess2D(crowdFF_model, pedestrianTransitionFunction);
	calAddElementaryProcess2D(crowdFF_model, conflictResolutionFunction);
	calAddElementaryProcess2D(crowdFF_model, updateDynamicFloorField);
	calRunAddSteeringFunc2D(crowdFF_simulation, clearConflicts);
	

	// simulation run
	calRunAddInitFunc2D(crowdFF_simulation, crowdFFInit);
	calRunInitSimulation2D(crowdFF_simulation);


	// Initialize the viewer
	calglInitViewer("crowdFF viewer", 1.0f, 400, 400, 40, 40, CAL_TRUE, 1);
	// model1 definition
	drawModel = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "pedestrians", crowdFF_model, crowdFF_simulation);
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
	
	drawModelConflicts = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "conflicts", crowdFF_model, crowdFF_simulation);
	calglAdd2Db (drawModelConflicts, NULL, &Q.conflicts, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_STATIC);
	calglAdd2Db (drawModelConflicts, Q.conflicts, &Q.conflicts, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_GREEN_SCALE, CALGL_DATA_TYPE_STATIC);
	
	drawModelFloorField = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "floor field", crowdFF_model, crowdFF_simulation);
	calglAdd2Dr (drawModelFloorField, NULL, &Q.floorField, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_STATIC);
	calglAdd2Dr (drawModelFloorField, Q.floorField, &Q.floorField, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_GREEN_SCALE, CALGL_DATA_TYPE_STATIC);
	
	drawModelDynamicFF = calglDefDrawModel2D (CALGL_DRAW_MODE_FLAT, "dynamicFF", crowdFF_model, crowdFF_simulation);
	calglAdd2Di (drawModelDynamicFF, NULL, &Q.dynamicFF, CALGL_TYPE_INFO_VERTEX_DATA, CALGL_TYPE_INFO_USE_NO_COLOR, CALGL_DATA_TYPE_STATIC);
	calglAdd2Di (drawModelDynamicFF, Q.dynamicFF, &Q.dynamicFF, CALGL_TYPE_INFO_COLOR_DATA, CALGL_TYPE_INFO_USE_GREEN_SCALE, CALGL_DATA_TYPE_STATIC);
	


	calglSetDisplayStep(1);

	calglMainLoop2D(argc, argv);

	

	return 0;
}

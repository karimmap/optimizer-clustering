#include "./limits.h"
#include "ortools/base/logging.h"
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <inttypes.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"


#include "problem.pb.h"
//#include "clusteringData.h"

using namespace std;
namespace operations_research {

void printProblem(problem::Problem& problem) {
  cout << "number of vehicles is " << problem.vehicles_size() << endl;
  cout << " number of services is " << problem.services_size() << endl;
  cout << "**************************Duration of each service "
          "is********************************"
       << endl;

  // for (auto service : problem.services())
  //   cout << service.duration() << " | " << endl;

  for (auto matrix : problem.matrices())
    cout << matrix.ShortDebugString() << endl;

  for (auto vehicle : problem.vehicles())
    cout << vehicle.DebugString() << endl;

  for (auto service : problem.services())
    cout << service.ShortDebugString() << endl;

  cout << "******************************************* vehicles capacity ****************************"
       << endl;
  for (auto vehicle : problem.vehicles()) {
    for (auto capacity : vehicle.capacities()) {
      cout << " V[" << vehicle.id() << "]=" << capacity.limit() << " | ";
    }
    cout << endl;
  }
  for( auto matrix: problem.matrices()){
    cout << " time.size() "<<matrix.time_size()<<endl;
  }
  cout << " ********************************** a vector of quantities for each sevice *************************************************"
       << endl;
  for (auto service : problem.services()) {
    for (int i = 0; i < service.quantities_size(); ++i){
      cout << " D[" << i << "]=" << service.quantities(i) << " | ";
    cout << endl;
  }
  }
}

void compatibleVehicle(problem::Problem& problem) {
  for (int i = 0; i < problem.services_size(); ++i) {
    auto service = problem.mutable_services(i);
    for (int k = 0; k < problem.vehicles_size(); ++k) {
      bool compatible = true;
      auto vehicle = problem.vehicles(k);
      for (int u = 0; u < service->quantities_size() && compatible; ++u) {
            auto quantity = service->quantities(u);
            auto capacity = vehicle.capacities(u);
            if (quantity > capacity.limit()) {
                compatible = false;
            }
       }
      if(service->skills_size() - vehicle.skills_size() > 0){
        compatible = false;
      }
        for(int s = 0; s < service->skills_size() && compatible; ++s) {
          for(int sv = 0; sv < vehicle.skills_size() && compatible; ++sv){
              if (service->skills(s) != vehicle.skills(s)){
                compatible = false;
              }
           }
        }

      if (compatible) {
        service->add_compatibale_vehicle_indices(k);
      }
    }
  }
}

void run() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  const string& filename = "instance-14.bin";
  problem::Problem problem;
  {
    fstream input(filename, ios::in | ios::binary);
    if (!problem.ParseFromIstream(&input)) {
      cout << "Failed to parse pbf." << endl;
    }
  }

  //printProblem(problem);

  compatibleVehicle(problem);

    for(auto service : problem.services()) {
     // for( int k = 0; k < service.compatibale_vehicle_indices_size(); ++k){
        cout << " "<< service.DebugString() <<endl;
     // }
    }

  // Create the linear solver with the GLOP backend.
  // MPSolver solver("simple_lp_program", MPSolver::GLOP_LINEAR_PROGRAMMING);
  MPSolver solver("simple_mip_program", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
  // Create the variables y.


  MPVariable*** const Y = new MPVariable**[problem.vehicles_size()];
  for (int k = 0; k < problem.vehicles_size(); ++k) {
    Y[k] = new MPVariable*[problem.services_size()];
    for (int j = 0; j < problem.services_size(); ++j) {
      stringstream ss;
      ss << "Y[" << k << "]"
         << "[" << j << "]";
      Y[k][j] = solver.MakeIntVar(0.0, 1.0, ss.str());
    }
  }
  MPVariable** const D = new MPVariable*[problem.vehicles_size()];
  for (int i = 0; i < problem.vehicles_size(); ++i) {
    stringstream ss;
    ss << "D[" << i << "]";
    D[i] = solver.MakeNumVar(0.0, INFINITY, ss.str());
  }

  MPVariable* const Dmax = solver.MakeNumVar(0.0, INFINITY, "Dmax");
  //---- Creation of the objective function-----//
  cout << "creating the objective function" << endl;

  MPObjective* const objective = solver.MutableObjective();
  objective->SetCoefficient(Dmax, 1);
  for (int k = 0; k < problem.vehicles_size(); k++) {
    objective->SetCoefficient(D[k], 1);
  }

  objective->SetMinimization();

  //----Creation of the constraints----//
  cout << "creating the constraint" << endl;

  //-----first constraint-----//

  for (int k = 0; k < problem.vehicles_size(); k++) {
    stringstream ss;
    ss << "";
    MPConstraint* const ct = solver.MakeRowConstraint(-INFINITY, 0.0, ss.str());
    ct->SetCoefficient(D[k], 1);
    ct->SetCoefficient(Dmax, -1);
  }

  // -----Second constraint ------//
  for(int i=0;i<problem.services_size();++i){
    for(int j=0;j<problem.services_size();++j){
      for(int k=0;k<problem.vehicles_size();++k){
        for(int m = 0; m < problem.matrices_size(); ++m ){
        stringstream ss;
        ss << "ToLin" <<"("<<i<<","<<j<<","<<k<<")";
        MPConstraint* const ct = solver.MakeRowConstraint(-problem.matrices(m).time(problem.services(i).matrix_index()),INFINITY,ss.str());
        ct->SetCoefficient(D[k],1);
        ct->SetCoefficient(Y[k][i],-problem.matrices(m).time(problem.services(i).matrix_index()));
        ct->SetCoefficient(Y[k][j],-problem.matrices(m).time(problem.services(i).matrix_index()));
       }
     }
   }
  }

  // assignment constraint
  for( int i = 0; i < problem.services_size(); ++i){
    auto service = problem.services(i);
    stringstream ss;
    ss << "assign("<<i<<")";
    MPConstraint* const ct = solver.MakeRowConstraint(1.0, 1.0, ss.str());
    for( int k =0; k < service.compatibale_vehicle_indices_size(); ++k  ){
           ct->SetCoefficient(Y[service.compatibale_vehicle_indices(k)][i],1);
    }
  }


  // capacity constraint for each vehicle
  for(int k = 0; k < problem.vehicles_size(); ++k){
    for(int u = 0; u < problem.vehicles(k).capacities_size(); ++u){
      stringstream ss;
      ss << "Cap(" <<k<< "," <<u<<")";
      if(problem.vehicles()[k].capacities()[u].limit() >= 0){
       MPConstraint* const ct = solver.MakeRowConstraint(0.0, problem.vehicles(k).capacities(u).limit(), ss.str());
         for(int i=0;i<problem.services_size();++i){
           ct->SetCoefficient(Y[k][i],problem.services(i).quantities(u));
         }
       }
    }
    cout << endl;
  }


  //-------respect avail time of vehicle--------------//

  for (int k = 0; k < problem.vehicles_size(); ++k) {
    stringstream ss;
    ss << "T(" << k << ")";
    MPConstraint* const ct = solver.MakeRowConstraint(0.0, problem.vehicles(k).duration(), ss.str());

    for (int i = 0;i < problem.services_size(); i++) {
      ct->SetCoefficient(Y[k][i], problem.services(i).duration());
    }
  }

  cout << "--> Configuring the solver" << endl;
  solver.set_time_limit(600); //< sets the time limit (in seconds)
  solver.SetNumThreads(1);    //< limits the solver to single thread usage
  string* lpFilePath = new string("model.lp");
  solver.ExportModelAsLpFormat(true, lpFilePath);

  // --- Solver launch ---
  cout << "--> Running the solver" << endl;
  const MPSolver::ResultStatus resultStatus = solver.Solve();

  // --- Solver results retrieval ---
  cout << "--> Retrieving solver results " << endl;
  if (resultStatus == MPSolver::OPTIMAL || (resultStatus == MPSolver::FEASIBLE)) {
    // the solver has computed the optimal solution or a feasible solution (when the time
    // limit is reached before proving optimality)
    cout << "Succes! (Status: " << resultStatus << ")"
         << endl; //< prints the solver status (see the documentation)
    cout << "Runtime : " << solver.wall_time() << " milliseconds" << endl;

    cout << "--> Printing results " << endl;
    cout << "Objective value = " << objective->Value()
         << endl; //<gets the value of the objective function for the best computed
                  //solution (optimal if no time limit)
    for (int k = 0; k < problem.vehicles_size(); ++k) {
      cout << "- Diameter " << k << " has " << D[k]->solution_value() << endl;
      for (int i = 0; i < problem.services_size(); ++i) {
        cout << "\t> service " << i << " is assigned to vehicle \t" << k << "-->"
             << Y[k][i]->solution_value() << endl;
      }
    }
    cout << "- Diameter Max   " << Dmax->solution_value() << endl;
  } else {
    // the model is infeasible (maybe wrong) or the solver has reached the time limit
    // without finding a feasible solution
    cerr << "Fail! (Status: " << resultStatus << ")"
         << endl; //< see status page in the documentation
  }

  delete lpFilePath;
  delete[] D;
  for (int i = 0; i < problem.vehicles_size(); ++i) {
    delete[] Y[i];
  }
  delete[] Y;
}


} // namespace operations_research

int main(/*int argc, char** argv*/) {
  operations_research::run();
  return EXIT_SUCCESS;
}
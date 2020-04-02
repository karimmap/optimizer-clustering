
#include "ortools/base/logging.h"
#include "ortools/linear_solver/linear_solver.h"
#include "ortools/linear_solver/linear_solver.pb.h"


#include <src/cpp/flann/flann.hpp>

#include "concave.h"

#include "problem.pb.h"

#include "ext-lib/json.h"

using json = nlohmann::json;

using namespace std;
namespace operations_research {

void printProblem(problem::Problem& problem) {
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
      auto vehicle_add = problem.mutable_vehicles(k);
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
        service->add_compatible_vehicle_indices(k);
      }
       vehicle_add -> set_id(k);
    }
    service -> set_id(i);
  }
}
 void readColor( vector < string > &color, string const & filePath ){
   cout << "<read> Read the data " << endl;
   ifstream f(filePath.c_str());
   if (!f.good()){
     cerr << "<read> Error: Could not open file " << filePath << "." << endl;
     exit(EXIT_FAILURE);
   }
   color.resize(20);
   for( int i = 0; i < 20; ++i){
     f >> color[i];
   }
   f.close();
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
        cout << " "<< service.DebugString() <<endl;
    }

  // MPSolver solver("simple_lp_program", MPSolver::GLOP_LINEAR_PROGRAMMING);
  MPSolver solver("simple_mip_program", MPSolver::CBC_MIXED_INTEGER_PROGRAMMING);
  const double infinity = solver.infinity(); // to use if a variable is not bounded or to create inequality constraints


 // Create the variables y.
  MPVariable*** const Y = new MPVariable**[problem.vehicles_size()];
  for (auto vehicle: problem.vehicles()) {
    Y[vehicle.id()] = new MPVariable*[problem.services_size()];
    for (auto service: problem.services()) {
      stringstream ss;
      ss << "Y[" << vehicle.id() << "]"
         << "[" << service.id() << "]";
      Y[vehicle.id()][service.id()] = solver.MakeIntVar(0.0, 1.0, ss.str());
    }
  }
  // Create the variables D.
  MPVariable** const D = new MPVariable*[problem.vehicles_size()];
  for (auto vehicle: problem.vehicles()) {
    stringstream ss;
    ss << "D[" << vehicle.id() << "]";
    D[vehicle.id()] = solver.MakeNumVar(0.0, infinity, ss.str());
  }
  // Create the variables Dmax
  MPVariable* const Dmax = solver.MakeNumVar(0.0, infinity, "Dmax");
  //---- Creation of the objective function-----//
  cout << "creating the objective function" << endl;

  MPObjective* const objective = solver.MutableObjective();
  objective->SetCoefficient(Dmax, 1);
  // for(auto vehicle: problem.vehicles()){
  //   objective->SetCoefficient(D[vehicle.id()], 1);
  // }

  objective->SetMinimization();

  //----Creation of the constraints----//
  cout << "creating the constraint" << endl;

  //-----first constraint-----//

  for (auto vehicle:problem.vehicles()) {
    stringstream ss;
    ss << "";
    MPConstraint* const ct = solver.MakeRowConstraint(-infinity, 0.0, ss.str());
    ct->SetCoefficient(D[vehicle.id()], 1);
    ct->SetCoefficient(Dmax, -1);
  }

 // -----Second constraint ------//
for(auto k : problem.vehicles()){
  auto matrix = problem.matrices(k.matrix_index()).distance();
  for(auto i : problem.services()){
    for(auto j : problem.services()){
      auto d_ijk = matrix.at(i.matrix_index() * sqrt(matrix.size()) + j.matrix_index());
      stringstream ss;
      ss << "ToLin" << "(" << i.id() << "," << j.id() << "," << k.id() << ")";
      MPConstraint* const ct = solver.MakeRowConstraint(-d_ijk, INFINITY, ss.str());
      ct->SetCoefficient(D[k.id()], 1);
      ct->SetCoefficient(Y[k.id()][i.id()], -d_ijk);
      ct->SetCoefficient(Y[k.id()][j.id()], -d_ijk);
    }
  }
}

  // assignment constraint
  for( int i = 0; i < problem.services_size(); ++i){
    auto service = problem.services(i);
    stringstream ss;
    ss << "assign("<<i<<")";
    MPConstraint* const ct = solver.MakeRowConstraint(1.0, 1.0, ss.str());
    for( int k =0; k < service.compatible_vehicle_indices_size(); ++k  ){
       ct->SetCoefficient(Y[service.compatible_vehicle_indices(k)][i],1);
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
  }

  //-------respect avail time of vehicle--------------//
  for(auto vehicle: problem.vehicles()) {
    stringstream ss;
    ss << "T(" << vehicle.id() << ")";
    MPConstraint* const ct = solver.MakeRowConstraint(0.0, vehicle.duration(), ss.str());
    for(auto service: problem.services()) {
       ct->SetCoefficient(Y[vehicle.id()][service.id()], service.duration());
    }
  }

  cout << "--> Configuring the solver" << endl;
  solver.set_time_limit(600); //< sets the time limit (in seconds)
  solver.SetNumThreads(1);    //< limits the solver to single thread usage

  string* modelInLpFormat = new string();
  solver.ExportModelAsLpFormat(false,modelInLpFormat);

  //print the model in a file
  ofstream writer ("model.lp");
  if(writer.is_open()) {
    writer << (*modelInLpFormat);
    writer.close();
  }
  else cout << "Unable to open the file and write the model";
  delete modelInLpFormat;

  // --- Solver launch ---
  cout << "--> Running the solver" << endl;
  const MPSolver::ResultStatus resultStatus = solver.Solve();

  // --- Solver results retrieval ---
  cout << "--> Retrieving solver results" << endl;
  if (resultStatus == MPSolver::OPTIMAL || (resultStatus == MPSolver::FEASIBLE)) {
    // the solver has computed the optimal solution or a feasible solution (when the time limit is reached before proving optimality)
    cout << "Succes! (Status: " << resultStatus << ")" << endl; //< prints the solver status (see the documentation)
    cout << "Runtime : " << solver.wall_time() << " milliseconds" << endl;
    cout << "--> Printing results " << endl;
    cout << "Objective value = " << objective->Value() << endl; //<gets the value of the objective function for the best computed
     //solution (optimal if no time limit)
    for (auto vehicle: problem.vehicles()) {
      cout << "- Diameter " << vehicle.id() << " has " << D[vehicle.id()]->solution_value() << endl;
      for (auto service: problem.services()) {
        cout << "\t> service " << service.id() << " is assigned to vehicle \t" << vehicle.id() << "-->"
             << Y[vehicle.id()][service.id()]->solution_value() << endl;
      }
    }
    cout << "- Diameter Max   " << Dmax->solution_value() << endl;
  }
  else {
    // the model is infeasible (maybe wrong) or the solver has reached the time limit without finding a feasible solution
    cerr << "Fail! (Status: " << resultStatus << ")" << endl; //< see status page in the documentation
  }

  //output solution
  vector < string > color;
  readColor(color,"color.txt");

  std::ofstream o("pretty.json");
  std::vector<json> v;

  for(int k = 0; k < problem.vehicles_size(); ++k){
    json t;
    PointVector points;
    for( int i = 0; i < problem.services_size(); ++i){
      for( int q = 0; q < problem.services(i).quantities_size(); ++q){
          if( Y[k][i]->solution_value() >= 1e-4){
            t["geometry"] = { {"type", "Point"},
            {"coordinates",
              { problem.services(i).location().longitude(),
              problem.services(i).location().latitude()}
              }
            };
            Point p;
            p.x = double(problem.services(i).location().longitude());
            p.y = double(problem.services(i).location().latitude());
            AddPoint(points, p);

            string lon_lat =  to_string(problem.services(i).location().latitude()) + "," + to_string(problem.services(i).location().longitude());
            string lat_lon =  to_string(problem.services(i).location().longitude()) + "," + to_string(problem.services(i).location().latitude());
            t["properties"] = { {"color",color[k]},{"stroke",color[k]},{"marker-size","large"},
            {"marker-color",color[k]},{"stroke-width", 10},{"name", "test_model"},
            {"lat_lon", lat_lon},
            {"lon_lat", lon_lat},
            {"duration",problem.vehicles(k).duration()},{"kg",55},
            {"qte",problem.services(i).quantities(q)},{"v-id",problem.vehicles(k).name()}, {"s-id",problem.services(i).name()} };
            t["type"] = "Feature";
            v.push_back(t);
          }
      }
    }

    bool iterate = true;
    int startpoint = 0;
    RemoveDuplicates(points);
    IdentifyPoints(points);
    startpoint = std::min(std::max(startpoint, 3), (int)points.size() - 1);
    PointVector hull = ConcaveHull(points, (size_t)startpoint, iterate);
    cout << endl;
    typedef double T;
    typedef std::array<T, 2> point_type;
    std::vector<point_type> Hull;
    for (int i = 0; i < hull.size(); i++){
      array<T, 2> m = {hull[i].x,hull[i].y};
      Hull.push_back(m);
    }
    int n = hull.size();
    if (Hull[0][0] != Hull[n - 1][0] && Hull[0][1] != Hull[n - 1][1] ){
      array<T, 2> m = {hull[0].x,hull[0].y};
      Hull.push_back(m);
    }
    t["geometry"] = { {"type", "Polygon"},
                     {"coordinates",{Hull} }
                    };
    string lon_lat =  to_string(problem.vehicles(k).start_location().latitude()) + "," + to_string(problem.vehicles(k).start_location().longitude());
    string lat_lon =  to_string(problem.vehicles(k).start_location().longitude()) + "," + to_string(problem.vehicles(k).start_location().latitude());
    t["properties"] = {{"color",color[k]},{"fill",color[k]},{"marker-size","small"},
              {"marker-color",color[k]},{"stroke-width", 10},{"name", problem.vehicles(k).name()},
              {"duration",problem.vehicles(k).duration()},{"v-id",problem.vehicles(k).name()},{"kg",55},
              {"qte",11},{"lat_lon", lat_lon},{"lon_lat", lon_lat} };
    v.push_back(t);

  }

  json s = { {"type","FeatureCollection"},{"features" , v} };
  o << std::setw(4) << s << std::endl;


  delete[] D;
  for (int i = 0; i < problem.vehicles_size(); ++i) {
    delete[] Y[i];
  }
  delete[] Y;
}

} // namespace operations_research

int main(int argc, char* argv[]) {
  operations_research::run();
  return EXIT_SUCCESS;

}
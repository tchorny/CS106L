#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <iterator>
#include "SimpleGraph.h"

using namespace std;

const double kPi = 3.14159265358979323;
const double kRepel = 0.001;
const double kAttract = 0.001;

void Welcome();
string GetLine(string);
void LoadGraph(SimpleGraph&, string);
void PopulateGraph(SimpleGraph&, ifstream&);
void ArrangeNodesOnCycle(SimpleGraph&, size_t);
vector<Node> ComputeRepelForce(const SimpleGraph&);
vector<Node> ComputeAttractForce(const SimpleGraph&);
void MoveGraph(SimpleGraph&, const vector<Node>&);
size_t ExtractSecondsFromString(const string&);

// Main method
int main() {
    Welcome();

    while (true) {

        string filename = GetLine("Please, choose graph filename from 'res' or enter 'quit'");
        if (filename == "quit") break;

        string timing = GetLine("Please, choose number of seconds for algorithm to run");
        const size_t seconds = ExtractSecondsFromString(timing);

        SimpleGraph graph;
        LoadGraph(graph, filename);
        InitGraphVisualizer(graph);
        DrawGraph(graph);

        time_t startTime = time(NULL);

        while (difftime(time(NULL), startTime) <= seconds) {

            MoveGraph(graph, ComputeRepelForce(graph));
            MoveGraph(graph, ComputeAttractForce(graph));
            DrawGraph(graph);

        }

    }

    return 0;
}

/* Prints a message to the console welcoming the user and
 * describing the program. */
void Welcome() {
    cout << "Welcome to CS106L GraphViz!" << endl;
    cout << "This program uses a force-directed graph layout algorithm" << endl;
    cout << "to render sleek, snazzy pictures of various graphs." << endl;
    cout << endl;
}

void LoadGraph(SimpleGraph& graph, string filename){

    ifstream input(filename.c_str());

    if (input) PopulateGraph(graph, input);
    else (cerr << "Couldn't open the file " << filename << endl);
}

string GetLine(string message = "") {
    if (!message.empty()) cout << message << endl;
    string filename;
    getline(cin, filename);
    return filename;
}

void PopulateGraph(SimpleGraph& graph, ifstream& input) {
    size_t nodes, index1, index2;

    if (input >> nodes) ArrangeNodesOnCycle(graph, nodes);

    while (input >> index1 >> index2) {
        Edge edge = {index1, index2};
        graph.edges.push_back(edge);
    }
}

void ArrangeNodesOnCycle(SimpleGraph& graph, size_t nodes) {
    for (size_t i = 0; i < nodes; ++i) {
        Node node = {cos(2*kPi*i/nodes), sin(2*kPi*i/nodes)};
        graph.nodes.push_back(node);
    }
}

vector<Node> ComputeRepelForce(const SimpleGraph& graph) {
    vector<Node> repelForce(graph.nodes.size(), {0, 0});
    for (auto iter1 = graph.nodes.begin(); iter1 != graph.nodes.end(); ++iter1) {
        for (auto iter2 = iter1 + 1; iter2 != graph.nodes.end(); ++iter2) {
            double force = kRepel / sqrt(pow(iter2->x - iter1->x, 2) + pow(iter2->y - iter1->y, 2));
            double theta = atan2(iter2->y - iter1->y, iter2->x - iter1->x);
            repelForce.at(iter1 - graph.nodes.begin()).x -= force * cos(theta);
            repelForce.at(iter1 - graph.nodes.begin()).y -= force * sin(theta);
            repelForce.at(iter2 - graph.nodes.begin()).x += force * cos(theta);
            repelForce.at(iter2 - graph.nodes.begin()).y += force * sin(theta);
        }

    }
    return repelForce;
}

vector<Node> ComputeAttractForce(const SimpleGraph& graph) {
    vector<Node> attractForce(graph.nodes.size(), {0, 0});
    for (auto iter = graph.edges.begin(); iter != graph.edges.end(); ++iter) {
        size_t start = iter->start;
        size_t end = iter->end;
        double x_0 = graph.nodes.at(start).x;
        double y_0 = graph.nodes.at(start).y;
        double x_1 = graph.nodes.at(end).x;
        double y_1 = graph.nodes.at(end).y;
        double force = kAttract * (pow(x_1 - x_0, 2) + pow(y_1 - y_0, 2));
        double theta = atan2(y_1 - y_0, x_1 - x_0);
        attractForce.at(start).x += force * cos(theta);
        attractForce.at(start).y += force * sin(theta);
        attractForce.at(end).x -= force * cos(theta);
        attractForce.at(end).y -= force * sin(theta);
    }
    return attractForce;
}

void MoveGraph(SimpleGraph& graph, const vector<Node>& displacement) {
    if(graph.nodes.size() == displacement.size()) {
        for (auto iter = graph.nodes.begin(); iter != graph.nodes.end(); ++iter) {
            iter->x += displacement.at(distance(graph.nodes.begin(), iter)).x;
            iter->y += displacement.at(distance(graph.nodes.begin(), iter)).y;
        }
    }
    else cerr << "Displacement vector doesn't match Graph nodes vector" << endl;
}

size_t ExtractSecondsFromString(const string& timing) {
    istringstream input (timing);
    size_t seconds;
    input >> seconds;
    return seconds;
}

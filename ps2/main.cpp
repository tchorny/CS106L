
#include <iostream>
#include <unordered_set>
#include <vector>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include "wikiscraper.h"


using std::cout;            using std::endl;
using std::string;          using std::vector;
using std::priority_queue;  using std::unordered_map;
using std::unordered_set;   using std::find;

void printVector(const vector<string>& vec) {
    for(auto item : vec) cout << item << " ";
    cout << endl;
}

/*
 * This is the function you will be implementing. It takes
 * two string representing the names of a start_page and
 * end_page and is supposed to return a ladder, represented
 * as a vector<string>, of links that can be followed from
 * start_page to get to the end_page.
 *
 * For the purposes of this algorithm, the "name" of a Wikipedia
 * page is what shows at the end of the URL when you visit that page
 * in your web browser. For ex. the name of the Stanford University
 * Wikipedia page is "Stanford_University" since the URL that shows
 * in your browser when you visit this page is:
 *
 *       https://en.wikipedia.org/wiki/Stanford_University
 */
vector<string> findWikiLadder(const string& start_page, const string& end_page) {

    WikiScraper scraper;

    auto target_set = scraper.getLinkSet(end_page);

    auto cmpFn = [&scraper, &target_set](vector<string> left, vector<string> right) {
        size_t leftPriority = 0, rightPriority = 0;
        auto leftSet = scraper.getLinkSet(left[left.size() - 1]);
        auto rightSet = scraper.getLinkSet(right[right.size() - 1]);
        for(auto item : target_set){
            if (leftSet.find(item) != leftSet.end()) ++leftPriority;
            if (rightSet.find(item) != rightSet.end()) ++rightPriority;
        }
        return leftPriority < rightPriority;
    };
    priority_queue<vector<string>, vector< vector<string> >, decltype(cmpFn)> pq(cmpFn);

    vector<string> vec{start_page};
    pq.push(vec);

    while(!pq.empty()) {
        vector<string> candidate(pq.top());
        pq.pop();
        printVector(candidate);
        auto candidate_set = scraper.getLinkSet(candidate[candidate.size() - 1]);
        if (candidate_set.find(end_page) != candidate_set.end()) {
            candidate.push_back(end_page);
            return candidate;
        }
        for (auto item : candidate_set) {
            if (find(candidate.begin(), candidate.end(), item) == candidate.end()) {
                vector<string> newCandidate(candidate);
                newCandidate.push_back(item);
                pq.push(newCandidate);
            }
        }
    }
    return {};
}

int main() {

    time_t startTime = time(NULL);

    auto ladder = findWikiLadder("Milkshake", "Gene");
    cout << endl;

    if(ladder.empty()) {
        cout << "No ladder found!" << endl;
    } else {
        cout << "Ladder found:" << endl;
        cout << "\t";

        printVector(ladder);
    }

    cout << "Time: " << difftime(time(NULL), startTime) << endl;

    return 0;
}





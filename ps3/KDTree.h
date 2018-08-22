/**
 * File: KDTree.h
 * Author: Aleksandr Chornyi (tchorny@gmail.com)
 * ------------------------
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 */

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>
#include <unordered_map>
#include <queue>
#include <utility>
#include <vector>

// "using namespace" in a header file is conventionally frowned upon, but I'm
// including it here so that you may use things like size_t without having to
// type std::size_t every time.
using namespace std;

template <size_t N, typename ElemType>
class KDTree {
public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();
    
    // Destructor: ~KDTree()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTree.
    ~KDTree();
    
    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Deep-copies the contents of another KDTree into this one.
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;
    
    // size_t size() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree and whether the tree is
    // empty.
    size_t size() const;
    bool empty() const;
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>& pt) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);
    
    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    
    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, const size_t k) const;

private:

    struct Node {
        Point<N> point;
        ElemType image;
        Node* left;
        Node* right;
    };

    Node* root;
    size_t numNodes;

    void clearSubtree(Node* node);

    //Helper functions to handle 'contains', '[]', 'at', 'insert'
    //methods. Returns nullptr if not fount.
    Node** findNode(const Point<N>& pt);
    Node const * const * findNode(const Point<N> &pt) const;

    void kNNValueHelper(BoundedPQueue<const Node*>&, const Point<N>& key,
                        const Node* const curr, const size_t pivotCoord) const;
    void copyOther(const KDTree& other);
    void copySubTree(const Node* const from, Node** to);
};

/** KDTree class implementation details */

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() : root(nullptr), numNodes(0) {
    if(N == 0) throw std::out_of_range("Zero-dimensional point");
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
    copyOther(rhs);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs){
    if(this != &rhs) {
        clearSubtree(root);
        copyOther(rhs);
    }
    return *this;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {clearSubtree(root);}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::clearSubtree(Node *node) {
    if(node != nullptr) {
        if(node->left != nullptr) clearSubtree(node->left);
        if(node->right != nullptr) clearSubtree(node->right);
        delete node;
    }
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::copyOther(const KDTree& other) {
    root = nullptr;
    numNodes = other.numNodes;
    copySubTree(other.root, &root);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::copySubTree(const Node* const from, Node** to) {
    if(from == nullptr) return;
    *to = new Node {from->point, from->image, nullptr, nullptr};
    copySubTree(from->left, &((*to)->left));
    copySubTree(from->right, &((*to)->right));
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {return N;}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {return numNodes;}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {return size() == 0;}

//Pointer to pointer helps to work with insertion to nullptr subtree
template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    Node** curr = findNode(pt);
    if(*curr == nullptr) {
        *curr = new Node {pt, value, nullptr, nullptr};
        ++numNodes;
    }
    else (*curr)->image = value;
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::Node const * const *
KDTree<N, ElemType>::findNode(const Point<N>& pt) const{
    Node const * const * curr =  &root;
    size_t pivotCoord = 0;
    while((*curr != nullptr) && ((*curr)->point != pt)) {
        if(pt[pivotCoord % N] < (*curr)->point[pivotCoord % N]) curr = &((*curr)->left);
        else curr = &((*curr)->right);
        ++pivotCoord;
    }
    return curr;
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::Node** KDTree<N, ElemType>::findNode(const Point<N>& pt) {
    Node** curr =  &root;
    size_t pivotCoord = 0;
    while((*curr != nullptr) && ((*curr)->point != pt)) {
        if(pt[pivotCoord % N] < (*curr)->point[pivotCoord % N]) curr = &((*curr)->left);
        else curr = &((*curr)->right);
        ++pivotCoord;
    }
    return curr;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {return *findNode(pt) != nullptr;}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
    Node** curr = findNode(pt);
    if(*curr == nullptr) {
        *curr = new Node {pt, ElemType(), nullptr, nullptr};
        ++numNodes;
    }
    ElemType& value = (*curr)->image;
    return value;
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    Node const * const * curr = findNode(pt);
    if(*curr == nullptr) throw std::out_of_range("Point is not in the tree or zero-dimension");
    return (*curr)->image;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
    Node** curr = findNode(pt);
    if(*curr == nullptr) throw std::out_of_range("Point is not in the tree or zero-dimension");
    return (*curr)->image;
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, const size_t k) const {

    if(k == 0) throw std::out_of_range("0-nearest neighbors search");
    BoundedPQueue<const Node*> bpq(k);
    kNNValueHelper(bpq, key, root, 0);

    if(bpq.empty()) throw std::out_of_range("Empty KDTree kNN search");

    unordered_map<ElemType, size_t> frequencyCounter;
    while(!bpq.empty()) ++frequencyCounter[bpq.dequeueMin()->image];

    auto cmp = [](pair<ElemType, size_t> p1, pair<ElemType, size_t> p2)
    { return p1.second < p2.second; };

    priority_queue<pair<ElemType, size_t>, vector< pair<ElemType, size_t> >,
            decltype(cmp)> imageCounterPQueue(cmp);

    for(auto item : frequencyCounter) imageCounterPQueue.push(item);

    return imageCounterPQueue.top().first;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::kNNValueHelper(BoundedPQueue<const Node*>& bpq, const Point<N>& key,
                                         const Node* const curr, const size_t pivotCoord) const {
    if(curr == nullptr) return;
    double distance = Distance(curr->point, key);
    if(distance < bpq.worst() || bpq.size() < bpq.maxSize()) bpq.enqueue(curr, distance);

    const double distanceToPivot = key[pivotCoord % N] - curr->point[pivotCoord % N];

    (distanceToPivot < 0) ? kNNValueHelper(bpq, key, curr->left, pivotCoord + 1) :
                            kNNValueHelper(bpq, key, curr->right, pivotCoord + 1);

    if(bpq.size() < bpq.maxSize() || fabs(distanceToPivot) < bpq.worst()) {

        (distanceToPivot < 0) ? kNNValueHelper(bpq, key, curr->right, pivotCoord + 1) :
                                kNNValueHelper(bpq, key, curr->left, pivotCoord + 1);
    }
}


#endif // KDTREE_INCLUDED

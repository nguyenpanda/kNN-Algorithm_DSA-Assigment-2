#include "main.hpp"
#include "Dataset.hpp"

/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */

struct kDTreeNode {
    vector<int> data;
    int label;
    kDTreeNode* left;
    kDTreeNode* right;

    kDTreeNode(vector<int> data, int label = 0, kDTreeNode* left = nullptr, kDTreeNode* right = nullptr) {
        this->data = data;
        this->label = label;
        this->left = left;
        this->right = right;
    }

    friend ostream& operator<<(ostream& os, const kDTreeNode& node) {
        os << "(";
        for (int i = 0; i < node.data.size(); i++) {
            os << node.data[i];
            if (i != node.data.size() - 1) {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }

    string change() {
        string s = "";
        s += "(";
        for (int i = 0; i < data.size(); i++) {
            s += to_string(data[i]);
            if (i == data.size() - 1) s += ")";
            else s += ", ";
        }
        return s;
    }
};

class kDTree {
private:
    int k;
    kDTreeNode* root;
    int count;

public:
    kDTree(int k = 2);

    kDTreeNode* copy(kDTreeNode* temp);

    void free(kDTreeNode* temp);

    ~kDTree();

    const kDTree& operator=(const kDTree& other);

    kDTree(const kDTree& other);

    void recInorderTraversal(kDTreeNode* temp, string& content) const;

    void recPreorderTraversal(kDTreeNode* temp, string& content) const;

    void recPostorderTraversal(kDTreeNode* temp, string& content) const;

    int recHeight(kDTreeNode* temp) const;

    int recLeaf(kDTreeNode* temp) const;

    kDTreeNode* newNode(const vector<int>& point, int label = 0);

    kDTreeNode* recInsert(kDTreeNode* temp, const vector<int>& point, int level);

    kDTreeNode* recRemove(kDTreeNode* temp, const vector<int>& point, int level);

    kDTreeNode* findMinNode(kDTreeNode* temp, int alpha, int level);

    kDTreeNode* compareThreeNode(kDTreeNode* x, kDTreeNode* left_x, kDTreeNode* right_x, int alpha);

    bool recSearch(kDTreeNode* temp, const vector<int>& point, int level);

    void merge(vector<vector<int>>& pointList, vector<int>& label, int start, int mid, int end, int u);

    void mergeSort(vector<vector<int>>& pointList, vector<int>& label, int start, int end, int u);

    kDTreeNode* recBuildTree(vector<vector<int>>& pointList, vector<int>& label, int level);

    double distance(const vector<int>& vt1, const vector<int>& vt2);

    void recNearestNeighbour(kDTreeNode* temp, const vector<int>& target, kDTreeNode*& best, int level, double& R_res);

    void sortBB(vector<kDTreeNode*>& bestList, const vector<int>& target);

    void recKNearestNeighbour(const vector<int>& target,
                              kDTreeNode* temp,
                              vector<kDTreeNode*>& bestList, int k_element,
                              int level);


    void inorderTraversal() const;

    void preorderTraversal() const;

    void postorderTraversal() const;

    int height() const;

    int nodeCount() const;

    int leafCount() const;

    void insert(const vector<int>& point);

    void remove(const vector<int>& point);

    bool search(const vector<int>& point);

    void buildTree(const vector<vector<int>>& pointList, vector<int> label = vector<int>());

    void nearestNeighbour(const vector<int>& target, kDTreeNode*& best);

    void kNearestNeighbour(const vector<int>& target, int k, vector<kDTreeNode*>& bestList);
};

class kNN {
private:
    int k;
    kDTree* treeKD;
    Dataset* X_train;
    Dataset* y_train;

public:
    kNN(int k = 5);

    ~kNN();

    void fit(Dataset& X_train, Dataset& y_train);

    Dataset predict(Dataset& X_test);

    double score(const Dataset& y_test, const Dataset& y_pred);
};

// Please add more or modify as needd

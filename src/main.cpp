#include "main.hpp"
#include "kDTree.hpp"

// tree.insert({5, 6}) syntax not support for macOS clang

void tc1() {
    kDTree tree(2);
    tree.insert(vector<int>({5, 6}));
    tree.insert(vector<int>({2, 2}));
    tree.insert(vector<int>({7, 3}));
    tree.insert(vector<int>({2, 8}));
    tree.insert(vector<int>({8, 7}));
    tree.insert(vector<int>({8, 1}));
    tree.insert(vector<int>({9, 4}));
    tree.insert(vector<int>({3, 5}));
    tree.inorderTraversal();
}

void tc2() {
    kDTree tree(2);
    tree.insert(vector<int>({5, 6}));
    tree.insert(vector<int>({2, 2}));
    tree.insert(vector<int>({7, 3}));
    tree.insert(vector<int>({2, 8}));
    tree.insert(vector<int>({8, 7}));
    tree.insert(vector<int>({8, 1}));
    tree.insert(vector<int>({9, 4}));
    tree.insert(vector<int>({3, 5}));
    tree.insert(vector<int>({9, 2}));
    tree.preorderTraversal();
}

//void tc3() {
//    kDTree tree(2);
//    vector<vector<int> > pointList = {{5, 6}, {2, 2}, {7, 3}, {2, 8}, {8, 7}, {8, 1}, {9, 4}, {3, 5}};
//    tree.buildTree(pointList);
//    tree.preorderTraversal();
//}
//
//void tc4() {
//    kDTree tree(2);
//    vector<vector<int> > pointList = {{5, 6}, {2, 2}, {7, 3}, {2, 8}, {8, 7}, {8, 1}, {9, 4}, {3, 5}};
//    tree.buildTree(pointList);
//    cout << tree.search({2, 2}) << endl;
//    cout << tree.search({9, 3}) << endl;
//}

void tc5() {
    kDTree tree(2);
    tree.insert(vector<int>({5, 6}));
    tree.insert(vector<int>({2, 2}));
    tree.insert(vector<int>({7, 3}));
    tree.insert(vector<int>({2, 8}));
    tree.insert(vector<int>({8, 1}));
    tree.insert(vector<int>({3, 5}));
    tree.insert(vector<int>({9, 2}));
    tree.remove(vector<int>({5, 6}));
    tree.preorderTraversal();
}

void tc6() {
    kDTree tree(2);
    tree.insert(vector<int>({5, 6}));
    tree.insert(vector<int>({2, 2}));
    tree.insert(vector<int>({7, 3}));
    tree.insert(vector<int>({2, 8}));
    tree.insert(vector<int>({3, 5}));
    tree.insert(vector<int>({8, 2}));
    tree.insert(vector<int>({8, 7}));
    tree.insert(vector<int>({9, 2}));
    tree.insert(vector<int>({9, 5}));
    kDTreeNode* best = nullptr;
    tree.nearestNeighbour({9, 3}, best);
    cout << "Nearest neighbour of (9, 3) is " << *best << endl;
}

void tc7() {
    kDTree tree(2);
    tree.insert(vector<int>({5, 6}));
    tree.insert(vector<int>({2, 2}));
    tree.insert(vector<int>({7, 3}));
    tree.insert(vector<int>({2, 8}));
    tree.insert(vector<int>({3, 5}));
    tree.insert(vector<int>({8, 2}));
    tree.insert(vector<int>({8, 7}));
    tree.insert(vector<int>({9, 2}));
    tree.insert(vector<int>({9, 5}));
    vector<kDTreeNode*> bestList;
    tree.kNearestNeighbour({9, 3}, 5, bestList);
    cout << "5 Nearest neighbour of (9, 3) are: ";
    for (auto node: bestList) {
        cout << *node << " ";
    }
}

void tc8() {
    Dataset dataset;
    dataset.loadFromCSV("mnist.csv");
    int nRows, nCols;

    kNN knn;
    Dataset X_train, X_test, y_train, y_test;
    Dataset feature = dataset.extract(0, -1, 1, -1);
    Dataset label = dataset.extract(0, -1, 0, 0);

    train_test_split(feature, label, 0.2, X_train, X_test, y_train, y_test);
    knn.fit(X_train, y_train);
    Dataset y_pred = knn.predict(X_test);

    std::cout << "y_pred" << endl;
    y_pred.printHead(10, 10);
    std::cout << endl;
    std::cout << "y_test" << endl;
    y_test.printHead(10, 10);
    std::cout << endl;
}

void tc9() {
    Dataset dataset;
    dataset.loadFromCSV("mnist.csv");
    int nRows, nCols;

    kNN knn;
    Dataset X_train, X_test, y_train, y_test;
    Dataset feature = dataset.extract(0, -1, 1, -1);
    Dataset label = dataset.extract(0, -1, 0, 0);

    train_test_split(feature, label, 0.2, X_train, X_test, y_train, y_test);
    knn.fit(X_train, y_train);
    Dataset y_pred = knn.predict(X_test);
    double accuracy = knn.score(y_test, y_pred);
    std::cout << "Accuracy: " << accuracy << endl;
}

int main(int argc, const char* argv[]) {
    tc1();

    return 0;
}
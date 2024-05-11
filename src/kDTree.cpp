#include "kDTree.hpp"

/* TODO: You can implement methods, functions that support your data structures here.
 * */

///////////////////// class kDTree /////////////////////
kDTree::kDTree(int k) {
    this->k = k;
    this->root = nullptr;
    this->count = 0;
}

kDTreeNode* kDTree::copy(kDTreeNode* temp) {
    if (temp == nullptr) return nullptr;
    kDTreeNode* res = new kDTreeNode(temp->data);
    res->label = temp->label;
    res->left = copy(temp->left);
    res->right = copy(temp->right);
    return res;
}

void kDTree::free(kDTreeNode* temp) {
    if (temp == nullptr) return;
    free(temp->left);
    free(temp->right);
    delete temp;
}

kDTree::~kDTree() {
    free(this->root);
    this->k = 0;
    this->count = 0;
}

const kDTree& kDTree::operator=(const kDTree& other) {
    if (this == &other) {
        return *this;
    }
    this->k = other.k;
    this->count = other.count;
    if (other.root != nullptr) {
        this->root = this->copy(other.root);
    } else {
        this->root = nullptr;
    }
    return *this;
}

kDTree::kDTree(const kDTree& other) {
    this->k = other.k;
    this->root = this->copy(other.root);
    this->count = other.count;
}

void kDTree::recInorderTraversal(kDTreeNode* temp, string& content) const {
    if (temp == nullptr) {
        return;
    }
    recInorderTraversal(temp->left, content);
    content += temp->change();
    content += " ";
    recInorderTraversal(temp->right, content);
}

void kDTree::inorderTraversal() const {
    string content = "";
    this->recInorderTraversal(root, content);
    if (!content.empty()) {
        content.pop_back();
    }
    cout << content;
}

void kDTree::recPreorderTraversal(kDTreeNode* temp, string& content) const {
    if (temp == nullptr) {
        return;
    }
    content += temp->change();
    content += " ";
    recPreorderTraversal(temp->left, content);
    recPreorderTraversal(temp->right, content);
}

void kDTree::preorderTraversal() const {
    string content = "";
    this->recPreorderTraversal(root, content);
    if (!content.empty()) {
        content.pop_back();
    }
    cout << content;
}

void kDTree::recPostorderTraversal(kDTreeNode* temp, string& content) const {
    if (temp == nullptr) {
        return;
    }
    recPostorderTraversal(temp->left, content);
    recPostorderTraversal(temp->right, content);
    content += temp->change();
    content += " ";
}

void kDTree::postorderTraversal() const {
    string content = "";
    this->recPostorderTraversal(root, content);
    if (!content.empty()) {
        content.pop_back();
    }
    cout << content;
}

int kDTree::height() const {
    return this->recHeight(root);
}

int kDTree::recHeight(kDTreeNode* temp) const {
    if (temp != nullptr) {
        int leftHeight = recHeight(temp->left);
        int rightHeight = recHeight(temp->right);
        return 1 + max(leftHeight, rightHeight);
    }
    return 0;
}

int kDTree::nodeCount() const {
    return this->count;
}

int kDTree::recLeaf(kDTreeNode* temp) const {
    if (temp == nullptr) {
        return 0;
    }
    if (temp->left == nullptr && temp->right == nullptr) {
        return 1;
    }
    int leftLeaf = recLeaf(temp->left);
    int rightLeaf = recLeaf(temp->right);
    return leftLeaf + rightLeaf;
}

int kDTree::leafCount() const {
    return this->recLeaf(root);
}

kDTreeNode* kDTree::newNode(const vector<int>& point, int label) {
    kDTreeNode* newNode = new kDTreeNode(point, label);
    return newNode;
}

kDTreeNode* kDTree::recInsert(kDTreeNode* temp, const vector<int>& point, int level) {
    if (temp == nullptr) {
        return newNode(point);
    }
    int alpha = level % k;
    if (point[alpha] < temp->data[alpha]) {
        temp->left = recInsert(temp->left, point, level + 1);
    } else {
        temp->right = recInsert(temp->right, point, level + 1);
    }
    return temp;
}

void kDTree::insert(const vector<int>& point) {
    this->root = recInsert(root, point, 0);
    this->count++;
}

kDTreeNode* kDTree::compareThreeNode(kDTreeNode* x, kDTreeNode* left_x, kDTreeNode* right_x, int alpha) {
    kDTreeNode* min_node = x;
    if (left_x && left_x->data[alpha] < min_node->data[alpha])
        min_node = left_x;
    if (right_x && right_x->data[alpha] < min_node->data[alpha])
        min_node = right_x;
    return min_node;
}

kDTreeNode* kDTree::findMinNode(kDTreeNode* temp, int alpha, int level) {
    if (temp == nullptr) {
        return nullptr;
    }
    int currentAlpha = level % k;
    if (currentAlpha == alpha) {
        if (temp->left == nullptr) {
            return temp;
        }
        return findMinNode(temp->left, alpha, level + 1);
    }
    return compareThreeNode(temp,
                            findMinNode(temp->left, alpha, level + 1),
                            findMinNode(temp->right, alpha, level + 1),
                            alpha);
}

kDTreeNode* kDTree::recRemove(kDTreeNode* temp, const vector<int>& point, int level) {
    if (temp == nullptr) {
        return nullptr;
    }

    int alpha = level % k;
    if (temp->data == point) {
        if (temp->right != nullptr) {
            kDTreeNode* min = findMinNode(temp->right, alpha, level + 1);
            temp->data = min->data;
            temp->right = recRemove(temp->right, min->data, level + 1);
        } else if (temp->left != nullptr) {
            kDTreeNode* min = findMinNode(temp->left, alpha, level + 1);
            temp->data = min->data;
            temp->right = temp->left;
            temp->left = nullptr;
            temp->right = recRemove(temp->right, min->data, level + 1);
        } else {
            delete temp;
            if (this->count > 0) {
                this->count--;
            }
            return nullptr;
        }
    } else if (point[alpha] < temp->data[alpha]) {
        temp->left = recRemove(temp->left, point, level + 1);
    } else {
        temp->right = recRemove(temp->right, point, level + 1);
    }
    return temp;
}

void kDTree::remove(const vector<int>& point) {
    this->root = recRemove(this->root, point, 0);
}

bool kDTree::recSearch(kDTreeNode* temp, const vector<int>& point, int level) {
    if (!temp || point.empty()) {
        return false;
    }

    if (temp->data == point) {
        return true;
    }

    int alpha = level % k;
    if (point[alpha] < temp->data[alpha]) {
        return recSearch(temp->left, point, level + 1);
    } else {
        return recSearch(temp->right, point, level + 1);
    }
}

bool kDTree::search(const vector<int>& point) {
    return recSearch(this->root, point, 0);
}

void kDTree::merge(vector<vector<int>>& pointList, vector<int>& label, int start, int mid, int end, int u) {
    int n1 = mid - start + 1;
    int n2 = end - mid;

    vector<vector<int>> leftPoint(n1), rightPoint(n2);
    vector<int> leftLabel(n1), rightLabel(n2);

    for (int i = 0; i < n1; i++) {
        leftPoint[i] = pointList[start + i];
        leftLabel[i] = label[start + i];
    }

    for (int j = 0; j < n2; j++) {
        rightPoint[j] = pointList[mid + 1 + j];
        rightLabel[j] = label[mid + 1 + j];
    }


    int x = 0, y = 0, z = start;

    while (x < n1 && y < n2) {
        if (leftPoint[x][u] <= rightPoint[y][u]) {
            pointList[z] = leftPoint[x];
            label[z] = leftLabel[x];
            x++;
        } else {
            pointList[z] = rightPoint[y];
            label[z] = rightLabel[y];
            y++;
        }
        z++;
    }

    while (x < n1) {
        pointList[z] = leftPoint[x];
        label[z] = leftLabel[x];
        x++;
        z++;
    }

    while (y < n2) {
        pointList[z] = rightPoint[y];
        label[z] = rightLabel[y];
        y++;
        z++;
    }
}

void kDTree::mergeSort(vector<vector<int>>& pointList, vector<int>& label, int start, int end, int u) {
    if (start < end) {
        int mid = start + (end - start) / 2;
        mergeSort(pointList, label, start, mid, u);
        mergeSort(pointList, label, mid + 1, end, u);
        merge(pointList, label, start, mid, end, u);
    }
}

kDTreeNode* kDTree::recBuildTree(vector<vector<int>>& pointList, vector<int>& label, int level) {

    if (pointList.empty() || pointList.size() == 0) {
        return nullptr;
    }

    int alpha = level % k;
    mergeSort(pointList, label, 0, pointList.size() - 1, alpha);
    int startIndex = 0;
    int endIndex = pointList.size() - 1;
    int mid = startIndex + (endIndex - startIndex) / 2;
    while (mid > 0 && pointList[mid - 1][alpha] == pointList[mid][alpha]) {
        mid--;
    }
    vector<vector<int>> left(pointList.begin(), pointList.begin() + mid);
    vector<vector<int>> right(pointList.begin() + mid + 1, pointList.end());
    vector<int> leftLabel(label.begin(), label.begin() + mid);
    vector<int> rightLabel(label.begin() + mid + 1, label.end());
    kDTreeNode* Node = newNode(pointList[mid], label[mid]);
    Node->left = recBuildTree(left, leftLabel, level + 1);
    Node->right = recBuildTree(right, rightLabel, level + 1);

    return Node;
}

void kDTree::buildTree(const vector<vector<int>>& pointList, vector<int> label) {
    if (label.empty()) {
        int n = pointList.size();
        for (int i = 0; i < n; i++) {
            label.push_back(0);
        }
    }
    vector<vector<int>> tem(pointList.begin(), pointList.end());
    this->root = recBuildTree(tem, label, 0);
    this->count = pointList.size();
}

double kDTree::distance(const vector<int>& vt1, const vector<int>& vt2) {
    long int sum = 0;
    for (int i = 0; i < vt1.size(); i++) {
        sum += pow(vt1[i] - vt2[i], 2);
    }
    return double(sqrt(sum));
}

void
kDTree::recNearestNeighbour(kDTreeNode* temp, const vector<int>& target, kDTreeNode*& best, int level, double& R_res) {

    if (temp == nullptr) {
        return;
    }

    int alpha = level % this->k;

    if (temp->data[alpha] > target[alpha])
        recNearestNeighbour(temp->left, target, best, level + 1, R_res);
    else {
        recNearestNeighbour(temp->right, target, best, level + 1, R_res);
    }

    double R_tem = distance(temp->data, target);
    if (temp->left == nullptr && temp->right == nullptr) {
        if (R_tem < R_res) {
            best = temp;
            R_res = R_tem;
        }
    } else if (R_tem < R_res) {
        best = temp;
        R_res = R_tem;
    }

    double r = abs(temp->data[alpha] - target[alpha]);

    if (R_res >= r) {
        if (temp->data[alpha] <= target[alpha])
            recNearestNeighbour(temp->left, target, best, level + 1, R_res);
        else {
            recNearestNeighbour(temp->right, target, best, level + 1, R_res);
        }
    }
}

void kDTree::nearestNeighbour(const vector<int>& target, kDTreeNode*& best) {
    if (root == nullptr) return;
    double R_res = 10000.0;
    recNearestNeighbour(this->root, target, best, 0, R_res);
}

void kDTree::sortBB(vector<kDTreeNode*>& bestList, const vector<int>& target) {
    int n = bestList.size();
    for (int i = 0; i < n - 1; ++i) {
        bool swapped = false;
        for (int j = 0; j < n - i - 1; ++j) {
            double dist1 = distance(bestList[j]->data, target);
            double dist2 = distance(bestList[j + 1]->data, target);
            if (dist1 > dist2) {
                swap(bestList[j], bestList[j + 1]);
                swapped = true;
            }
        }
        // Nếu không có sự hoán đổi nào trong một lượt duyệt, tức là mảng đã được sắp xếp, dừng vòng lặp
        if (!swapped) {
            break;
        }
    }
}

void kDTree::recKNearestNeighbour(const vector<int>& target,
                                  kDTreeNode* temp, vector<kDTreeNode*>& bestList,
                                  int k_element,
                                  int level) {
    if (temp == nullptr) return;
    int alpha = level % this->k;
    if (target[alpha] >= temp->data[alpha]) {
        recKNearestNeighbour(target, temp->right, bestList, k_element, level + 1);
    } else {
        recKNearestNeighbour(target, temp->left, bestList, k_element, level + 1);
    }

    double R_current = distance(temp->data, target);
    if (bestList.size() < k_element) {
        bestList.push_back(temp);
    } else {
        double maxR = distance(bestList[0]->data, target);
        int indexmaxR = 0;
        int n = bestList.size();
        for (int i = 1; i < n; i++) {
            double kc = distance(bestList[i]->data, target);
            if (kc > maxR) {
                maxR = kc;
                indexmaxR = i;
            }
        }

        if (R_current < maxR) {
            bestList[indexmaxR] = temp;
        }
        sortBB(bestList, target);
    }

    double R = distance(bestList.back()->data, target);
    int r = abs(target[alpha] - temp->data[alpha]);
    if (R >= r) {
        if (target[alpha] >= temp->data[alpha]) {
            recKNearestNeighbour(target, temp->left, bestList, k_element, level + 1);
        } else {
            recKNearestNeighbour(target, temp->right, bestList, k_element, level + 1);
        }
    }
}

void kDTree::kNearestNeighbour(const vector<int>& target, int k, vector<kDTreeNode*>& bestList) {
    if (root == nullptr) return;
    recKNearestNeighbour(target, root, bestList, k, 0);
}

///////////////////// class kNN /////////////////////
kNN::kNN(int k) {
    this->k = k;
    this->treeKD = new kDTree(this->k);
    this->X_train = nullptr;
    this->y_train = nullptr;
}

void kNN::fit(Dataset& X_train, Dataset& y_train) {
    this->X_train = &X_train;
    this->y_train = &y_train;

    vector<vector<int>> pointList;
    for (auto tem: this->X_train->data) {
        vector<int> x;
        for (int i: tem) {
            x.push_back(i);
        }
        pointList.push_back(x);
    }

    vector<int> label;
    for (auto i: this->y_train->data) {
        int x = i.front();
        label.push_back(x);
    }

    if (this->treeKD != nullptr) {
        delete this->treeKD;
    }
    this->treeKD = new kDTree(X_train.data.size());
    this->treeKD->buildTree(pointList, label);
}

Dataset kNN::predict(Dataset& X_test) {
    Dataset y_pred;
    vector<vector<int>> pointList;
    for (auto tem: this->X_train->data) {
        vector<int> x;
        for (int i: tem) {
            x.push_back(i);
        }
        pointList.push_back(x);
    }

    for (auto point: X_test.data) {
        vector<int> tem(point.begin(), point.end());
        vector<kDTreeNode*> bestList;
        this->treeKD->kNearestNeighbour(tem, k, bestList);

        vector<int> decodeLabel;
        for (int i = 0; i < bestList.size(); i++) {
            decodeLabel.push_back(bestList[i]->label);
        }

        int arrLabel[10];
        for (int i = 0; i < 10; i++) {
            arrLabel[i] = 0;
        }
        for (int i = 0; i < decodeLabel.size(); i++) {
            arrLabel[decodeLabel[i]]++;
        }
        int res_label = 0;
        for (int i = 0; i < 10; i++) {
            if (arrLabel[i] > arrLabel[res_label]) {
                res_label = i;
            }
        }
        list<int> label;
        label.push_back(res_label);
        y_pred.data.push_back(label);
    }
    y_pred.columnName.push_back("label");
    return y_pred;
}

double kNN::score(const Dataset& y_test, const Dataset& y_pred) {
    vector<vector<int>> vt_test;
    for (auto it: y_test.data) {
        vector<int> tem;
        for (int x: it) {
            tem.push_back(x);
        }
        vt_test.push_back(tem);
    }

    vector<vector<int>> vt_pred;
    for (auto it: y_pred.data) {
        vector<int> tem;
        for (int x: it) {
            tem.push_back(x);
        }
        vt_pred.push_back(tem);
    }

    int cnt = 0;
    int n = vt_test.size();

    for (int i = 0; i < n; i++) {
        if (vt_test[i][0] == vt_pred[i][0]) cnt++;
    }
    return double(cnt) / double(n);
}

kNN::~kNN() {
    if (this->treeKD != nullptr) {
        delete this->treeKD;
        this->treeKD = nullptr;
    }
}
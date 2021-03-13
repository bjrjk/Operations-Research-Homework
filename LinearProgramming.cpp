#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>

using namespace std;

class LinearProgramming {
private:
	const double MAX_D = 1e10; //大M
	int n, m; // n个变量，m个方程
	int relaxedN; // 松弛后relaxedN个变量
	int manualN; // 松弛后再添加人工变量，有manualN个变量
	vector<vector<double>> A; // A矩阵
	vector<double> b; //b矩阵
	vector<char> equationType; //方程类型
	vector<int> baseVar; //可行基下标是否存在及其对应的行号
	vector<double> sigmaArray; //检验数σ数组
	vector<double> thetaArray; //θ数组
	int solutionType; //0:还未运算完毕,1:无可行解,2:无穷多最优解,3:唯一最优解,4:无界解
	int bringInVarIndex, bringOutVarIndex; //换入变量下标，换出变量下标
	long scanLong(const string& str, int& index) { //扫描一个长整型
		char* ptr;
		long result = strtol(&str.c_str()[index], &ptr, 0);
		int newIndex = ptr - str.c_str();
		if (index == newIndex) {
			result = 1;
		}
		index = newIndex;
		return result;
	}
	double scanDouble(const string& str, int& index) { //扫描一个浮点数
		char* ptr;
		double result = strtod(&str.c_str()[index], &ptr);
		int newIndex = ptr - str.c_str();
		if (index == newIndex) {
			result = 1;
		}
		index = newIndex;
		return result;
	}
	void scanComponent(const string& str, int line, int& i, bool& position) { //扫描一个项
		while (i < str.size() && isspace(str[i]))i++; //除去空格
		if (i >= str.size())return;
		if (str[i] == '=' || str[i] == '<' || str[i] == '>') { //获取大小比较号
			equationType[line] = str[i];
			position = true;
			i++;
			return;
		}
		double coefficient = scanDouble(str, i); //扫描系数，如果x前无系数，自动置为1
		if (i >= str.size() || str[i] != 'x') { //单个数值，不带未知量
			b[line] += (position ? coefficient : -coefficient);
			return;
		}
		i++;
		long varIndex = scanLong(str, i); //扫描未知变量下标
		A[line][varIndex] += (position ? -coefficient : +coefficient);
	}
	void scanEquation(const string& str, int line) { //扫描一个方程/不等式
		int i = 0;
		bool position = false; //在大小比较号的左边为false，在右面为true
		while (i < str.size()) {
			scanComponent(str, line, i, position);
		}
	}
	string preprocessEquation(const string& str) { //类似于-x1 +x1这种x前面什么都不带的项，在中间插入一个1
		string s = "";
		for (int i = 0; i < str.size(); i++) {
			s += str[i];
			if (i + 1 < str.size()) {
				if ((str[i] == '-' || str[i] == '+') && str[i + 1] == 'x') {
					s += "1";
				}
			}
		}
		return s;
	}
	void calculateSigma() { //计算检验数σ
		for (int i = 0; i < manualN; i++) {
			if (!baseVar[i]) { //对于非基变量，计算检验数
				sigmaArray[i] = A[0][i]; //Cj
				for (int j = 0; j < manualN; j++) {
					if (baseVar[j]) { //减去基变量的Ci×Aij
						sigmaArray[i] -= A[0][j] * A[baseVar[j]][i];
					}
				}
			}
			else { //对于基变量，检验数为0
				sigmaArray[i] = 0;
			}
		}
	}
	bool checkSolution() { //检查可行解情况
		bool flag = true;
		for (int i = 0; i < manualN; i++) {
			if (sigmaArray[i] > 0) {
				flag = false;
				break;
			}
		}
		if (flag) { //所有σj<=0
			for (int i = relaxedN; i < manualN; i++) { //基变量中有非零的人工变量？
				if (baseVar[i]) {
					solutionType = 1; //无可行解
					return true;
				}
			}
			for (int i = 0; i < manualN; i++) { //非基变量检验数为0？
				if (!baseVar[i]) {
					if (fabs(sigmaArray[i]) < 1e-6) { //绝对值做差比较
						solutionType = 2; //无穷多最优解
						return true;
					}
				}
			}
			solutionType = 3; //唯一最优解
			return true;
		}
		else {
			for (int i = 0; i < manualN; i++) {
				if (sigmaArray[i] > 0) {
					flag = true;
					for (int j = 1; j <= m; j++) {
						if (A[j][i] > 0) {
							flag = false;
							break;
						}
					}
					if (flag) {
						solutionType = 4; //无界解
						return true;
					}
				}
			}
			solutionType = 0;//未计算完毕
			return false;
		}
	}
	void calculateBringIOVarIndex() { //计算换入换出变量下标
		double maxSigmaV = sigmaArray[0];
		int maxSigmaIndex = 0;
		for (int i = 1; i < manualN; i++) {
			if (sigmaArray[i] > maxSigmaV) {
				maxSigmaIndex = i;
				maxSigmaV = sigmaArray[i];
			}
		}
		bringInVarIndex = maxSigmaIndex; //确定换入变量
		double minThetaV = MAX_D;
		int minThetaIndex = -1;
		for (int i = 0; i < manualN; i++) {
			if (baseVar[i]) { //基变量
				if (A[baseVar[i]][bringInVarIndex] > 0) { //Aik>0
					thetaArray[baseVar[i]] = b[baseVar[i]] / A[baseVar[i]][bringInVarIndex];
					if (thetaArray[baseVar[i]] < minThetaV) {
						minThetaIndex = i;
						minThetaV = thetaArray[baseVar[i]];
					}
				}
				else {
					thetaArray[baseVar[i]] = NAN;
				}
			}
		}
		bringOutVarIndex = minThetaIndex;
	}
	void updateMatrix() { //迭代运算：替换基变量，更新矩阵元素
		baseVar[bringInVarIndex] = baseVar[bringOutVarIndex]; //换入变量
		baseVar[bringOutVarIndex] = 0; //换出变量
		vector<vector<double>> tmpMatrix = A;
		int mainLine = baseVar[bringInVarIndex], mainCol = bringInVarIndex;
		double Ask = tmpMatrix[mainLine][mainCol];
		double bs = b[mainLine];
		for (int i = 1; i <= m; i++) { //m行方程
			//操作b矩阵
			if (i == mainLine)b[i] /= Ask;
			else b[i] -= bs / Ask * tmpMatrix[i][mainCol];
			//操作A矩阵
			for (int j = 0; j < manualN; j++) { //manualN个变量
				if (i == mainLine) {
					A[i][j] /= Ask;
				}
				else if (j == mainCol) {
					A[i][j] = 0;
				}
				else {
					A[i][j] -= tmpMatrix[mainLine][j] / Ask * tmpMatrix[i][mainCol];
				}
			}
		}
	}
public:
	LinearProgramming(int n, int m) :n(n), m(m) { //构造函数
		A.resize(m + 1); //第0行放待求函数
		for (int i = 0; i <= m; i++)
			A[i].resize(n);
		b.resize(m + 1);
		equationType.resize(m + 1);
		baseVar.resize(n);
		thetaArray.resize(m + 1);
		solutionType = 0;
	}
	void init(vector<string>& equationArray) { //初始化方程
		assert(equationArray.size() == m + 1);
		for (int i = 0; i <= m; i++) {
			scanEquation(preprocessEquation(equationArray[i]), i);
		}
	}
	void printOriginEquations() { //打印化简后的原始方程
		cout << "\n";
		cout << "Print Simplified Equations:\n";
		for (int i = 0; i <= m; i++) { // m个方程依次打印
			for (int j = 0; j < n; j++) { // n个变量依次打印
				cout << A[i][j] << "x" << j << "\t";
			}
			if (i != 0) // 待求函数不打印b值
				cout << equationType[i] << "\t" << b[i];
			else
				cout << "=\tz -> max";
			cout << "\n";
		}
		cout << "Number of variable: " << n << "\n";
		cout << "Number of equation: " << m << "\n";
		cout << "\n";
	}
	void relax() { //添加松弛变量
		relaxedN = n;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] != '=')relaxedN++;
		}
		for (int i = 0; i <= m; i++)
			A[i].resize(relaxedN);
		baseVar.resize(relaxedN);
		int cnt = n - 1;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] != '=') {
				cnt++;
				if (equationType[i] == '<') { //添加正松弛变量
					A[i][cnt] = 1;
					baseVar[cnt] = i; //添加到初始可行基，行号为i
				}
				else { //equationType[i] == '>' //添加负松弛变量
					A[i][cnt] = -1;
				}
			}
		}
	}
	void printRelaxedEquations() { //打印松弛后的方程
		cout << "Print Relaxed Equations:\n";
		for (int i = 0; i <= m; i++) { // m个方程依次打印
			for (int j = 0; j < relaxedN; j++) { // n个变量依次打印
				cout << A[i][j] << "x" << j << "\t";
			}
			if (i != 0) // 待求函数不打印b值
				cout << "=\t" << b[i];
			else
				cout << "=\tz -> max";
			cout << "\n";
		}
		cout << "Number of relaxed variable: " << relaxedN << "\n";
		cout << "Number of relaxed equation: " << m << "\n";
		cout << "\n";
	}
	void manual() { //添加人工变量
		manualN = relaxedN;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] == '=' || equationType[i] == '>')manualN++;
		}
		for (int i = 0; i <= m; i++)
			A[i].resize(manualN);
		baseVar.resize(manualN);
		for (int i = relaxedN; i < manualN; i++) {
			A[0][i] = -MAX_D; //设置人工变量系数为-M
		}
		int cnt = relaxedN - 1;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] == '=' || equationType[i] == '>') { //添加人工变量
				cnt++;
				A[i][cnt] = 1;
				baseVar[cnt] = i; //添加到初始可行基，行号为i
			}
		}
	}
	void printManualEquations() { //打印松弛、添加人工变量的方程
		cout << "Print Manual Equations:\n";
		for (int i = 0; i <= m; i++) { // m个方程依次打印
			for (int j = 0; j < manualN; j++) { // n个变量依次打印
				cout << A[i][j] << "x" << j << "\t";
			}
			if (i != 0) // 待求函数不打印b值
				cout << "=\t" << b[i];
			else
				cout << "=\tz -> max";
			cout << "\n";
		}
		cout << "Number of manual variable: " << manualN << "\n";
		cout << "Number of manual equation: " << m << "\n";
		cout << "Initial Feasible Base: ";
		for (int i = 0; i < manualN; i++) {
			if (baseVar[i]) {
				cout << i << " ";
			}
		}
		cout << "\n\n";
	}
	void printSolveProcedure(int step) {
		cout << "Print Solve Procedure " << step << ":\n";
		cout << "Cj\t\t\t";
		for (int i = 0; i < manualN; i++)cout << A[0][i] << "\t";
		cout << "\n";
		cout << "CB\tXB\tb\t";
		for (int i = 0; i < manualN; i++)cout << "x" << i << "\t";
		cout << "θ\t\n";
		for (int i = 0; i < manualN; i++) {
			if (baseVar[i]) {
				cout << A[0][i] << "\t";
				cout << "x" << i << "\t";
				cout << b[baseVar[i]] << "\t";
				for (int j = 0; j < manualN; j++) {
					if (bringOutVarIndex == i && bringInVarIndex == j)
						cout << "[";
					cout << setprecision(4) << A[baseVar[i]][j];
					if (bringOutVarIndex == i && bringInVarIndex == j)
						cout << "]";
					cout << "\t";
				}
				if (!isnan(thetaArray[baseVar[i]]))
					cout << thetaArray[baseVar[i]] << "\t";
				cout << "\n";
			}
		}
		cout << "Cj-Zj\t\t\t";
		for (int i = 0; i < manualN; i++) {
			cout << sigmaArray[i] << "\t";
		}
		cout << "\n\n";
	}
	void printResult() {
		switch (solutionType) {
			case 1: {
				cout << "无可行解\n";
				break;
			}
			case 2: {
				cout << "无穷多最优解，一个最优解如下：\n";
				printSingleAnsResult();
				break;
			}
			case 3: {
				cout << "唯一最优解：\n";
				printSingleAnsResult();
				break;
			}
			case 4: {
				cout << "无界解\n";
				break;
			}
		}
	}
	void printSingleAnsResult() {
		double zValue = 0;
		for (int i = 0; i < manualN; i++) {
			if (baseVar[i]) {
				cout << "x" << i << "=" << b[baseVar[i]] << " ";
				zValue += A[0][i] * b[baseVar[i]];
			}
			else {
				cout << "x" << i << "=0 ";
			}
		}
		cout << "\n";
		cout << "max(z)= " << zValue << "\n";
	}
	void solve() {
		sigmaArray.resize(manualN);
		for (int i = 1; solutionType == 0; i++) {
			bringInVarIndex = bringOutVarIndex = -1;
			for (int i = 0; i < thetaArray.size(); i++)thetaArray[i] = NAN;
			calculateSigma();
			if (checkSolution()) {
				printSolveProcedure(i);
				break;
			}
			calculateBringIOVarIndex();
			printSolveProcedure(i);
			updateMatrix();
		}
		printResult();
	}
};
int main() {
	cout << "用大M法单纯型法求解线性规划问题，x变量下标从0开始：\n";
	int n, m;
	cout << "Input n(number of variable), m(number of equation): ";
	cin >> n >> m;
	cin.get();
	LinearProgramming LP(n, m);
	vector<string> equationArray;
	string s;
	cout << "Input z = ";
	getline(cin, s);
	equationArray.push_back(s);
	for (int i = 1; i <= m; i++) {
		cout << "Input equation " << i << ": ";
		getline(cin, s);
		equationArray.push_back(s);
	}
	LP.init(equationArray);
	LP.printOriginEquations();
	LP.relax();
	LP.printRelaxedEquations();
	LP.manual();
	LP.printManualEquations();
	LP.solve();
}
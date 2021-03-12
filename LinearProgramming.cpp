#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cctype>
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
	vector<int> baseVar; //可行基下标
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
public:
	LinearProgramming(int n, int m) :n(n), m(m) { //构造函数
		A.resize(m + 1); //第0行放待求函数
		for (int i = 0; i <= m; i++)
			A[i].resize(n);
		b.resize(m + 1);
		equationType.resize(m + 1);
	}
	void init(vector<string>& equationArray) { //初始化方程
		assert(equationArray.size() == m + 1);
		for (int i = 0; i <= m; i++) {
			scanEquation(preprocessEquation(equationArray[i]), i);
		}
	}
	void printOriginEquations() { //打印化简后的原始方程
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
	}
	void relax() { //添加松弛变量
		relaxedN = n;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] != '=')relaxedN++;
		}
		for (int i = 0; i <= m; i++)
			A[i].resize(relaxedN);
		int cnt = n - 1;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] != '=') {
				cnt++;
				if (equationType[i] == '<') { //添加正松弛变量
					A[i][cnt] = 1;
					baseVar.push_back(cnt); //添加到初始可行基
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
	}
	void manual() { //添加人工变量
		manualN = relaxedN;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] == '=' || equationType[i] == '>')manualN++;
		}
		for (int i = 0; i <= m; i++)
			A[i].resize(manualN);
		for (int i = relaxedN; i < manualN; i++) {
			A[0][i] = -MAX_D; //设置人工变量系数为-M
		}
		int cnt = relaxedN - 1;
		for (int i = 1; i <= m; i++) {
			if (equationType[i] == '=' || equationType[i] == '>') { //添加人工变量
				cnt++;
				A[i][cnt] = 1;
				baseVar.push_back(cnt);
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
	}
};
int main() {
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
}
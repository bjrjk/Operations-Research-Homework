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
	int n, m; // n个变量，m个方程
	vector<vector<double>> A;
	vector<double> b;
	vector<char> equationType;
	long scanLong(const string& str, int& index) {
		char* ptr;
		long result = strtol(&str.c_str()[index], &ptr, 0);
		int newIndex = ptr - str.c_str();
		if (index == newIndex) {
			result = 1;
		}
		index = newIndex;
		return result;
	}
	double scanDouble(const string& str, int& index) {
		char* ptr;
		double result = strtod(&str.c_str()[index], &ptr);
		int newIndex = ptr - str.c_str();
		if (index == newIndex) {
			result = 1;
			//newIndex++;
		}
		index = newIndex;
		return result;
	}
	void scanComponent(const string& str, int line, int& i, bool& position) {
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
	void scanEquation(const string& str, int line) {
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
	LinearProgramming(int n, int m) :n(n), m(m) {
		A.resize(m + 1); //第0行放待求函数
		for (int i = 0; i <= m; i++)
			A[i].resize(n);
		b.resize(m + 1);
		equationType.resize(m + 1);
	}
	void init(vector<string>& equationArray) {
		assert(equationArray.size() == m + 1);
		for (int i = 0; i <= m; i++) {
			scanEquation(preprocessEquation(equationArray[i]), i);
		}
	}
	void printEquations() {
		cout << "Print Equations:\n";
		for (int i = 0; i <= m; i++) { // m个方程依次打印
			for (int j = 0; j < n; j++) { // n个变量依次打印
				cout << A[i][j] << "x" << j << "\t";
			}
			if (i != 0) // 待求函数不打印b值
				cout << "=\t" << b[i];
			else
				cout << "=\tz";
			cout << "\n";
		}
	}
};
int main() {
	int n, m;
	cin >> n >> m;
	cin.get();
	LinearProgramming LP(n, m);
	vector<string> equationArray;
	for (int i = 0; i <= m; i++) {
		string s;
		getline(cin, s);
		equationArray.emplace_back(s);
	}
	LP.init(equationArray);
	LP.printEquations();
}
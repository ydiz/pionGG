#pragma once

#include <Grid/Grid.h>

// read pion pion correlator
std::vector<double> read_pp(const std::string &filename) {
	std::vector<double> ret;
	std::ifstream f(filename);

	double tmp;
	while(f) {
		f >> tmp;
		if(f.eof()) break;
		ret.push_back(tmp);

		f >> tmp;
		if(f.peek() == ',') f.ignore();
	}
	f.close();
	return ret;
}

namespace Grid {

namespace QCD {

template<class T>
void writeScidac(T& field, const std::string &filename){ // because of writeScidacFieldRecord, field cannot be const
  emptyUserRecord record;
  ScidacWriter WR(field._grid->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
  WR.open(filename);
	WR.writeScidacFieldRecord(field, record);
  WR.close();
};

template<class T>
void readScidac(T& field, const std::string &filename){
  emptyUserRecord record;
  ScidacReader RD;
  RD.open(filename);
	RD.readScidacFieldRecord(field, record);
  RD.close();
};




}}

#ifndef _AUX_H
#define _AUX_H
#include <bitset>
#include <limits>
#include <tuple>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

namespace aux{
template<typename VectorType>
std::pair<double,double> map_geometry(const std::vector<std::tuple<size_t,VectorType,double>>& flippable, double& decayCoeff){
    double minDistance = std::numeric_limits<double>::max();
    double maxDistance = 0.0;
    for(const auto& atom1 : flippable){
        for(const auto& atom2 : flippable){
             auto d = VectorType::norm(std::get<1>(atom1) - std::get<1>(atom2));
            if (d > 1e-7 && d < minDistance) minDistance = d;
            if (d > maxDistance)             maxDistance = d;
        }
    }
    decayCoeff = 24*log10(2)*log(10)/(maxDistance-minDistance);
    return std::make_pair(minDistance,maxDistance);
}

template<typename VectorType>
struct System{
    std::array<VectorType,3> basis;
    std::vector<std::tuple<size_t,VectorType,double>> supercell; 
    std::vector<std::tuple<size_t,VectorType,double>> flippable; 
};
    
template<typename VectorType>
std::vector<std::tuple<unsigned,unsigned,double>> get_interactions(System<VectorType>& system, double decayCoeff, size_t interactionModel=1U){
#ifdef _VERBOSE
	switch(interactionModel){
	  case 0U:
	    std::cout<<"Model of interactions: "<<"exponential"<<std::endl;
	    break;
	  case 1U:
	    std::cout<<"Model of interactions: "<<"rational"<<std::endl;
	    break;
	  case 2U:
	    std::cout<<"Model of interactions: "<<"Bessel J0"<<std::endl;
	    break;
	  default: exit(-1); break;
	}
#endif
    std::vector<std::tuple<unsigned,unsigned,double>> d;
    int i, j, k;
    double J = 0.0;
    VectorType image;
    for (const auto& atom1 : system.flippable) {
        for (const auto& atom2 : system.supercell) {
            J = 0.0;
            for (i = -1; i < 2; ++i) for (j = -1; j < 2; ++j) for (k = -1; k < 2; ++k) {
                image = i*system.basis[0] + j*system.basis[1] + k*system.basis[2];
                if (std::fabs(std::get<2>(atom2)) > 1e-9 && std::get<0>(atom1) != std::get<0>(atom2))
                    switch (interactionModel) {
                        case 0:
                            J += -exp(-decayCoeff*VectorType::norm(std::get<1>(atom1)-std::get<1>(atom2)-image));
                            break;
                        case 1:
                            J += -pow(VectorType::norm(std::get<1>(atom1)-std::get<1>(atom2)-image),-decayCoeff);
                            break;
                        case 2:
                            J += -sin(-decayCoeff*VectorType::norm(std::get<1>(atom1)-std::get<1>(atom2)-image))
                                      /decayCoeff*VectorType::norm(std::get<1>(atom1)-std::get<1>(atom2)-image);
                            break;
                        default: exit(-1); break;
                    }
            }
            d.push_back(std::make_tuple(std::get<0>(atom1), std::get<0>(atom2), J));
        }
    }
    return d;
}

template<typename State,size_t SITESNUMBER>
void print_state(const State& state, size_t reference, std::bitset<SITESNUMBER> mask){
    for (unsigned i=0; i<reference; ++i){
		if(mask[i])  std::cout<<"\033[32m"<<state[i]<<"\033[39m";
		else         std::cout<<state[i];
    }
    std::cout<<"\033[1m"<<state[reference]<<"\033[0m";
    for (unsigned i=reference+1; i<SITESNUMBER; ++i){
		if(mask[i])  std::cout<<"\033[91m"<<state[i]<<"\033[39m";
		else         std::cout<<state[i];
    }
	    std::cout<<std::endl;
}

template<typename VectorType>
void read_basis(std::array<VectorType,3>& basis,
                char _basisFile[]){ 
    std::size_t pos,buff;
    std::string line;
    std::ifstream directions(_basisFile);
    if(directions.is_open()) {
        size_t id = 0;
        while (getline(directions,line)){
            buff = 0;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            basis[id++] = position;
        }
        directions.close();
    }
}

template<typename VectorType>
void read_supercell(std::vector<std::tuple<size_t,VectorType,double>>& supercell, char _supercellFile[]){ 
    std::size_t pos,buff;
    std::string line;
    std::ifstream inSupercell(_supercellFile);
    if(inSupercell.is_open()) {
        while (getline(inSupercell,line)){
            size_t idx = std::stoi(line,&pos);
            buff = pos;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            auto moment = std::stod(line.substr(buff));
            supercell.push_back(std::make_tuple(idx,position,moment));
        }
        inSupercell.close();
    }
}

template<typename VectorType>
void read_flippables(std::vector<std::tuple<size_t,VectorType,double>>& flippable, char _flippableFile[]){ 
    std::size_t pos,buff;
    std::string line;
    std::ifstream inFlippable(_flippableFile);
    if(inFlippable.is_open()) {
        while (getline(inFlippable,line)){
            size_t idx = std::stoi(line,&pos);
            buff = pos;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            auto moment = std::stod(line.substr(buff));
            flippable.push_back(std::make_tuple(idx,position,moment));
        }
        inFlippable.close();
    }
}
void greetings();
} //end of namespace aux
#endif

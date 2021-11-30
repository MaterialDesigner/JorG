#include "solver.h"

int solver(char _basis[],char _supercell[],char _flippable[], size_t reference, size_t unique_flips, size_t ansatz){
#ifndef _SITESNUMBER
#define _SITESNUMBER 64
#endif
constexpr size_t SITESNUMBER = _SITESNUMBER;
   
// *****************************************************************************************
// **************************************   READ   *****************************************
// *****************************************************************************************
    aux::System<typename ising::IsingModel<SITESNUMBER>::VectorType> system;
    aux::read_basis(system.basis,_basis);
    aux::read_supercell(system.supercell,_supercell);
    aux::read_flippables(system.flippable,_flippable);

#ifndef _QUIET
    aux::greetings();
#endif

// *****************************************************************************************
// *******************************  INITIALIZE ENGINE    ***********************************
// *****************************************************************************************
    std::random_device randomDevice{};
    std::mt19937 generator{randomDevice()};
    std::normal_distribution<> gauss{0.0,1.0};

    ising::IsingModel<SITESNUMBER> model;

    model.set_basis(system.basis);
    model.set_supercell(system.flippable);
    model.set_reference(reference);

// *****************************************************************************************
// ********************************   MAP INTERACTIONS   ***********************************
// *****************************************************************************************
    double decayCoeff;
    aux::map_geometry(system.flippable,decayCoeff);

    std::bitset<SITESNUMBER> mask(std::string(SITESNUMBER,'0'));

    for (const auto& atom : system.flippable) mask.set(std::get<0>(atom));
    mask.reset(reference);

// *****************************************************************************************
// ********************************   CHECK MODEL   ****************************************
// *****************************************************************************************

#ifdef _VERBOSE    
            std::cout<<"               ";
            for(int i=1; i<(SITESNUMBER-reference); ++i) std::cout<<" ";
            std::cout<<"| reference"<<std::endl;
            std::cout<<"Mask:          "<<mask<<std::endl;
#endif

#ifdef _LOWERLIMIT
    size_t lowerlimit = _LOWERLIMIT;
#else
    size_t lowerlimit = 0.15*mask.count(); // doesn't allow less than 15% flipps
#endif
#ifdef _UPPERLIMIT
    size_t upperlimit = _UPPERLIMIT;
#else
    size_t upperlimit = 0.85*mask.count(); // doesn't allow more than 85% flipps
#endif
    upperlimit = upperlimit < 2U         ? 2U : upperlimit;
    lowerlimit = lowerlimit > upperlimit ? 1U : lowerlimit;
    lowerlimit = lowerlimit < 1U         ? 1U : lowerlimit;

    std::unordered_set<std::bitset<SITESNUMBER>> solutions;
    size_t iteration = 0U;
    for(auto n = 4.0; n<1024.0; n*=2){
        model.reset();
        auto d = aux::get_interactions(system,n*decayCoeff,ansatz);
        model.add_interaction(d);

#ifdef _VERBOSE    
        std::cout<<"Limits are set to be: ["<<lowerlimit<<","<<upperlimit<<']'<<std::endl;
#endif

        for(unsigned m = 0; m<unique_flips; ++m){
            model.randomize_state();

            auto x = model.run(&mask);
	    if(x.count() < lowerlimit) continue;
	    if(x.count() > upperlimit) continue;
	    if(x.count() > 0.5*mask.count()) x.flip();

	    auto isUnique = solutions.insert(x);

        if(isUnique.second){
#ifndef _QUIET
          std::cout<<"("<<iteration<<")\t"<<n*decayCoeff<<"\t";
          aux::print_state(x,reference,mask);
#endif
          ++iteration;
        }
        else {
#ifndef _QUIET
          std::cout<<"-------------------------------------------------------------"<<std::endl;
          aux::print_state(x,reference,mask);
#endif
        }
	    if(solutions.size() >= unique_flips) break;
        }
    if(solutions.size() >= unique_flips) break;
    }

    std::ofstream ostrm("best.flips", std::ios::out);
    for(const auto& x : solutions){
        for(unsigned b=0U; b<SITESNUMBER; ++b){
          ostrm<<x[b]<<" ";
        }
        ostrm<<std::endl;
    }

    return 0;
}

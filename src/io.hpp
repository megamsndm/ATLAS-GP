#pragma once
#include "system.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <map>

inline void write_xyz(const std::string& path, const System& sys,
                      double energy, const std::string& comment="") {
    std::ofstream f(path);
    if (!f) { std::cerr<<"Cannot write: "<<path<<"\n"; return; }
    f<<sys.N<<"\n"<<std::fixed<<std::setprecision(10)
     <<"Energy="<<energy<<" Formula="<<sys.formula();
    if(!comment.empty()) f<<" "<<comment;
    f<<"\n";
    for(int i=0;i<sys.N;++i)
        f<<sys.species[sys.type_id[i]]
         <<std::setw(20)<<sys.px(i)
         <<std::setw(20)<<sys.py(i)
         <<std::setw(20)<<sys.pz(i)<<"\n";
}

inline void write_xyz_append(const std::string& path, const System& sys,
                              double energy, int index=0) {
    std::ofstream f(path,std::ios::app);
    if(!f)return;
    f<<sys.N<<"\n"<<std::fixed<<std::setprecision(10)
     <<"Isomer="<<index<<" Energy="<<energy
     <<" Formula="<<sys.formula()<<"\n";
    for(int i=0;i<sys.N;++i)
        f<<sys.species[sys.type_id[i]]
         <<std::setw(20)<<sys.px(i)
         <<std::setw(20)<<sys.py(i)
         <<std::setw(20)<<sys.pz(i)<<"\n";
}

inline void write_isomers(const std::string& dir,
                          const std::string& formula,
                          const std::vector<std::pair<double,System>>& isomers,
                          const std::string& pot_name) {
    std::filesystem::create_directories(dir);
    std::string traj=dir+"/"+formula+"_isomers.xyz";
    std::string best=dir+"/"+formula+"_global_min.xyz";
    std::ofstream clear(traj,std::ios::trunc); clear.close();
    for(int i=0;i<(int)isomers.size();++i)
        write_xyz_append(traj,isomers[i].second,isomers[i].first,i);
    if(!isomers.empty())
        write_xyz(best,isomers[0].second,isomers[0].first,"Global_min Pot="+pot_name);
    std::cout<<"[ATLAS] "<<isomers.size()<<" isomers → "<<dir<<"\n";
}

inline System read_xyz(const std::string& path) {
    std::ifstream f(path);
    if(!f) throw std::runtime_error("Cannot read: "+path);
    int N; f>>N;
    std::string comment; std::getline(f,comment); std::getline(f,comment);
    System sys; sys.N=N; sys.pos.resize(3*N); sys.type_id.resize(N);
    std::map<std::string,int> smap;
    for(int i=0;i<N;++i){
        std::string el; double x,y,z; f>>el>>x>>y>>z;
        if(smap.find(el)==smap.end()){smap[el]=(int)sys.species.size();sys.species.push_back(el);}
        sys.type_id[i]=smap[el];
        sys.pos[3*i]=x;sys.pos[3*i+1]=y;sys.pos[3*i+2]=z;
    }
    sys.charges.assign(N,0.0); sys.mol_id.assign(N,-1);
    return sys;
}

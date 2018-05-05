#ifndef POOL_ELEMENT_H
#define POOL_ELEMENT_H

#include <map>
#include <string>
#include <iostream>

namespace CPL {


template <class T >
class Pool: public std::map<std::string, T*> {
public:
    using std::map<std::string, T*>::map;
    virtual ~Pool(){std::cout << "POOL DESTRUCTOR" << std::endl; clearPool();};
    virtual void clearPool(){
        typename std::map<std::string, T*>::iterator f;
        //TODO: Make this work properly. Fails with OpenFOAM producing a SegFault.
        // Free the memmory
        // for (f = this->begin(); f != this->end(); f++)
        //     delete f->second;
        // this->clear();

    }
};

template <class T >
class PoolElement {
public:
    PoolElement(){};
    virtual ~PoolElement(){};
    PoolElement(std::string name): elem_name(name){};
    std::string elem_name;
    Pool<T>* pool;

    virtual void addToPool(Pool<T>* T_pool) {
        //TODO: Use exceptions
        if(!T_pool->insert(std::pair<std::string, T*>(elem_name, static_cast<T*>(this))).second) {
            std::string err_msg = "RUNTIME-ERROR: Dependency already exists.";
            std::cout << err_msg << std::endl;
        }
        pool = T_pool;
    }
};
}
#endif //POOL_ELEMENT_H

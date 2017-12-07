#ifndef POOL_ELEMENT_H
#define POOL_ELEMENT_H

#include <map>
#include <string>

namespace CPL {

template <class T >
class PoolElement {
public:
    PoolElement(){};
    PoolElement(std::string name): elem_name(name){};
    std::string elem_name;
    std::map<std::string, T*>* pool;

    void addToPool(std::map<std::string, T*>& T_pool) {
        //TODO: Use exceptions
        if(!T_pool.insert(std::pair<std::string, T*>(elem_name, static_cast<T*>(this))).second) {
            std::string err_msg = "Dependency already exists.";
        }
        pool = &T_pool;
    }
};
}
#endif //POOL_ELEMENT_H

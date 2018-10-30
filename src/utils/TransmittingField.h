#ifndef TRANSMITTING_FIELD_H 
#define TRANSMITTING_FIELD_H

#include <vector>
#include <valarray>
#include <map>

#include "PoolElement.h"
#include "CPL_field.h"
#include "CPL_ndArray.h"

namespace CPL{

class BufferPtr {
    public:
        BufferPtr(): buffOffset(0) {}
        BufferPtr(CPL::ndArray<double>& buff_ptr, int buff_offset) :
                  buffPtr(&buff_ptr), buffOffset(buff_offset) {}
        int buffOffset;
        CPL::ndArray<double>* buffPtr;
        double& operator () (int i0, int i1, int i2, int i3) {
            return (*buffPtr)(buffOffset+i0, i1, i2, i3);
        }
        const double operator () (int i0, int i1, int i2, int i3) const {
            return (*buffPtr)(buffOffset+i0, i1, i2, i3);
        }
};



class TransmittingField: public CPL::PortionField, public CPL::PoolElement<TransmittingField> {
public:
    TransmittingField() {}
    TransmittingField(std::string field_name, const CPL::PortionField& portion_field, const CPL::PortionField& field) : 
                  PoolElement(field_name), CPL::PortionField(portion_field), region(field) {}
    TransmittingField(std::string field_name) : PoolElement(field_name) {}
    virtual void setup()=0;
    virtual void update(){}
    // virtual void addToPool(std::map<std::string, TransmittingField*>* pool, bool overwrite=true);
    virtual void addToPool(CPL::Pool<TransmittingField>* pool, bool overwrite=true);
    void setBuffer(const BufferPtr& buff) {buffer = buff;}
    int data_size;
    CPL::PortionField region; 
    CPL::BufferPtr buffer;
    virtual ~TransmittingField(){};
protected:
    void iteratePortion_();
    virtual void iterateFunc_(const std::vector<int>& glob_cell,
                              const std::vector<int>& loc_cell, 
                              const std::valarray<double>& coord)=0;
    virtual void iterateFuncCustom_()=0;
 
};

class IncomingField: public CPL::TransmittingField {
public:
    using TransmittingField::TransmittingField;
    void unpack() {iteratePortion_();}
    virtual ~IncomingField(){};
protected:
    void iterateFunc_(const std::vector<int>& glob_cell,
                      const std::vector<int>& loc_cell, 
                      const std::valarray<double>& coord) 
                      {unpack_(glob_cell, loc_cell, coord);}
    void iterateFuncCustom_() {unpackCustom_();}
    virtual void unpack_(const std::vector<int>& glob_cell, 
                        const std::vector<int>& loc_cell,
                        const std::valarray<double>& coord) {}
    virtual void unpackCustom_() {}
};

class OutgoingField: public CPL::TransmittingField {
public:
    using TransmittingField::TransmittingField;
    void pack() {iteratePortion_();}
    virtual ~OutgoingField(){};
protected:
    void iterateFunc_(const std::vector<int>& glob_cell,
                      const std::vector<int>& loc_cell, 
                      const std::valarray<double>& coord) 
                      {pack_(glob_cell, loc_cell, coord);}
    void iterateFuncCustom_() {packCustom_();}
    virtual void pack_(const std::vector<int>& glob_cell, 
                       const std::vector<int>& loc_cell,
                       const std::valarray<double>& coord) {}
    virtual void packCustom_() {}
 
};


// typedef std::map<std::string, TransmittingField*> TransmittingFieldMap;
typedef CPL::Pool<TransmittingField> TransmittingFieldPoolT;

class TransmittingFieldPool : public CPL::Pool<TransmittingField> {
    public:
        // TransmittingFieldPool(): TransmittingFieldMap() {}
        TransmittingFieldPool(): CPL::Pool<TransmittingField>() {}
        TransmittingFieldPool(const CPL::PortionField& portion_region,
                              const CPL::PortionField& region) : 
                              // TransmittingFieldMap(), portionField(portion_region), 
                              CPL::Pool<TransmittingField>(), portionField(portion_region), 
                              field(region) {}
        void setupAll();
        void updateAll();
        void allocateBuffer(CPL::ndArray<double>& buffer) const;
        CPL::ndArray<double> buffer;
        CPL::PortionField portionField;
        CPL::PortionField field;
    protected:
};

class OutgoingFieldPool : public TransmittingFieldPool {
    public:
        using CPL::TransmittingFieldPool::TransmittingFieldPool;
        void packAll();
};

class IncomingFieldPool : public TransmittingFieldPool {
    public:
        using CPL::TransmittingFieldPool::TransmittingFieldPool;
        void unpackAll();
};


}

#endif //TRANSMITTING_FIELD_H

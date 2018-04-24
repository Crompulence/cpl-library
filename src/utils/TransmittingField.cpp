#include "TransmittingField.h"
#include "cpl.h"

void CPL::TransmittingFieldPool::setupAll() {
    CPL::Pool<TransmittingField>::iterator f;
    for (f = this->begin(); f != this->end(); f++)
        f->second->setup();
}


void CPL::TransmittingFieldPool::updateAll() {
    CPL::Pool<TransmittingField>::iterator f;
    for (f = this->begin(); f != this->end(); f++)
        f->second->update();
}

void CPL::IncomingFieldPool::unpackAll() {
    CPL::Pool<TransmittingField>::iterator f;
    for (f = this->begin(); f != this->end(); f++)
        static_cast<CPL::IncomingField*>(f->second)->unpack();
}

void CPL::OutgoingFieldPool::packAll() {
    CPL::Pool<TransmittingField>::iterator f;
    for (f = this->begin(); f != this->end(); f++)
        static_cast<CPL::OutgoingField*>(f->second)->pack();
}


// void CPL::TransmittingField::addToPool(std::map<std::string, TransmittingField*>* 
void CPL::TransmittingField::addToPool(CPL::Pool<TransmittingField>* 
                                       pool, bool overwrite) {
    // Add portion and region information to the transmitting field
    CPL::PortionField a, b;
    a =  (static_cast<TransmittingFieldPool*>(pool))->portionField;
    b = (static_cast<TransmittingFieldPool*>(pool))->field;
    if (overwrite) {
        static_cast<PortionField&>(*this) = a;
        region = b;
    }
    CPL::PoolElement<TransmittingField>::addToPool(pool);
}


void CPL::TransmittingFieldPool::allocateBuffer(CPL::ndArray<double>& buffer) const {
    CPL::PortionField send_portion, recv_portion;
    // for (f = this.begin(); f != this.end(); f++) {
    //     if this->second->
    //
    // Receive buffer allocation
    int buff_size = 0;
    CPL::TransmittingFieldPool::const_iterator f;
    for (f = this->begin(); f != this->end(); f++) {
        f->second->setBuffer(CPL::BufferPtr(buffer, buff_size));
        buff_size += f->second->data_size;
    }
    int buff_shape[4] = {buff_size, portionField.nCells[0], portionField.nCells[1], portionField.nCells[2]};
    buffer.resize (4, buff_shape);

}

void CPL::TransmittingField::iteratePortion_() {
    std::vector<int> portion = cellBounds;
	if (CPL::is_proc_inside(portion.data())) {
        // Pack/unpack custom called here
        iterateFuncCustom_();
        std::vector<int> glob_cell(3), loc_cell(3);
        std::valarray<double> cell_coord(3);
		for (int i = portion[0]; i <= portion[1]; i++)
        for (int j = portion[2]; j <= portion[3]; j++)
        for (int k = portion[4]; k <= portion[5]; k++) {
            glob_cell = std::vector<int>({i, j, k});
            loc_cell = getLocalCell(glob_cell);
			// CPL::map_glob2loc_cell(portion.data(), glob_cell.data(), 
            //                        loc_cell.data());
            CPL::map_cell2coord(glob_cell[0], glob_cell[1], glob_cell[2], &cell_coord[0]); 
            // Pack/unpack called here
            iterateFunc_(glob_cell, loc_cell, cell_coord);
        }
	}
}

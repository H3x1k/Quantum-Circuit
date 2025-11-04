// Index implementation

#include "Index.hpp"

#include <iostream>

Index::Index(size_t i) : i({ i }) {}

Index::Index(std::vector<size_t> mi) : i(std::move(mi)) {}

Index::Index(std::initializer_list<size_t> il) : i(il) {}

Index Index::range(size_t s, size_t e) {
	std::vector<size_t> mi;
	for (int i = s; i < e; i++)
		mi.push_back(i);
	return Index(mi);
}

size_t& Index::operator[](size_t idx) {
	if (idx >= i.size())
		throw std::out_of_range("Index out of bounds");
	return i[idx];
}

const size_t& Index::operator[](size_t idx) const {
	if (idx >= i.size())
		throw std::out_of_range("Index out of bounds");
	return i[idx];
}
// Index implementation

#include "Index.hpp"

Index::Index(size_t i) : i({ i }) {}

Index::Index(std::vector<size_t> mi) : i(mi) {}

Index Index::range(size_t s, size_t e) {
	std::vector<size_t> mi;
	for (int i = s; i < e; i++)
		mi.push_back(i);
	return Index(mi);
}
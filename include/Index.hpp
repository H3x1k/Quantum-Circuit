// Index header
#pragma once

#include <vector>

class Index {
public:
	std::vector<size_t> i;

	Index(size_t i);
	Index(std::vector<size_t> mi);
	Index(std::initializer_list<size_t> il);
	static Index range(size_t s, size_t e);
	size_t& operator[](size_t idx);
	const size_t& operator[](size_t idx) const;
	//static Index all();
};
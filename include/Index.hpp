// Index header

#include <vector>

class Index {
	std::vector<size_t> i;
public:
	Index(size_t i);
	Index(std::vector<size_t> mi);
	static Index range(size_t s, size_t e);
	//static Index all();
};
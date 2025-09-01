// Index header

#include <vector>

class Index {
	bool multiIndex;
	size_t i;
	std::vector<size_t> mi;
public:
	Index(size_t i);
	Index(std::vector<size_t> mi);
};
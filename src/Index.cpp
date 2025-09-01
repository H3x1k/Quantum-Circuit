// Index implementation

#include "Index.hpp"

Index::Index(size_t i) : i(i), multiIndex(false) {}

Index::Index(size_t mi) : mi(mi), multiIndex(true) {}
#ifndef TRIE_H_
#define TRIE_H_

#include <map>
#include <fstream>

struct Trie_node {
  std::map<int, Trie_node*> children_map;

  ~Trie_node() {
    for (auto map_pair: children_map)
      delete map_pair.second;
  }
};

template <typename Alcove_id>
class Trie {

  std::size_t size_of_subtree(Trie_node* node) {
    std::size_t result = 0;
    if (node->children_map.empty())
      result++;
    for (auto m_pair: node->children_map)
      result += size_of_subtree(m_pair.second);
    return result;
  }
  
  Trie_node root_;
  
public:
  Trie() {}

  void add(const Alcove_id& a_id) {
    Trie_node* curr_node = &root_;
    for (auto& c: a_id) {
      auto map_it = curr_node->children_map.find(c);
      if (map_it == curr_node->children_map.end()) {
        Trie_node* new_node = new Trie_node;
        curr_node->children_map.emplace(c, new_node);
        curr_node = new_node;
      }
      else
        curr_node = map_it->second;
    }
  }

  bool contains(const Alcove_id& a_id) {
    Trie_node* curr_node = &root_;
    for (auto& c: a_id) {
      auto map_it = curr_node->children_map.find(c);
      if (map_it == curr_node->children_map.end())
        return false;
      else
        curr_node = map_it->second;
    }
    return true;
  }

  std::size_t size() {
    return size_of_subtree(&root_);
  }
  
  Trie_node const& root() const {
    return root_;
  }
  
};

std::ostream& operator<<(std::ostream & os, const Trie_node& trie_node) {
  os << "[ ";
  if (!trie_node.children_map.empty()) {
    auto map_it = trie_node.children_map.begin();
    os << map_it->first << *map_it->second;
    ++map_it;
    for (; map_it != trie_node.children_map.end(); ++map_it) {
      os << ", " << map_it->first << *map_it->second;
    }
  }
  os << " ]";
  return os;
}

template<class Alcove_id>
std::ostream& operator<<(std::ostream & os, const Trie<Alcove_id>& trie) {
  os << trie.root();
  return os; 
}

#endif

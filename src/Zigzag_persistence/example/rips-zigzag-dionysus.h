#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>

#include <boost/range/adaptors.hpp>
namespace ba = boost::adaptors;

#include <dionysus/simplex.h>
#include <dionysus/distances.h>
#include <dionysus/rips.h>
namespace d = dionysus;

#include <dionysus/dlog/progress.h>

typedef         std::vector<float>                                      Point;
typedef         std::vector<Point>                                      PointContainer;

typedef         d::PairwiseDistances<PointContainer,
                                     d::L2Distance<Point>>              PairDistances;
typedef         PairDistances::DistanceType                             DistanceType;
typedef         PairDistances::IndexType                                Vertex;

typedef         d::Rips<PairDistances>                                  Generator;
typedef         Generator::Simplex                                      Simplex;
typedef         std::set<Simplex>                                       SimplexSet;

typedef         std::vector<Vertex>                                     VertexVector;
typedef         std::vector<DistanceType>                               EpsilonVector;
typedef         std::tuple<Vertex,Vertex>                               Edge;
typedef         std::vector<Edge>                                       EdgeVector;

inline PointContainer compute_points(unsigned int numberOfPoints, int seed = -1)
{
	PointContainer finalPoints;
	std::set<Point> points;
	std::random_device dev;
	std::mt19937 rng(dev());
	if (seed > -1) rng.seed(seed);
	std::uniform_real_distribution<float> dist(0,10);

	for (unsigned int i = 0; i < numberOfPoints; ++i){
		auto res = points.insert({dist(rng), dist(rng)});
		while(!res.second){
			res = points.insert({dist(rng), dist(rng)});
		}
		finalPoints.push_back(*res.first);
	}

	return finalPoints;
}

inline void compute_vertices_and_epsilons(const PairDistances& distances, 
								   VertexVector& vertices, 
								   EpsilonVector& epsilons) 
{
	DistanceType inf = std::numeric_limits<DistanceType>::infinity();
	EpsilonVector dist(distances.size(), inf);

	vertices.push_back(distances.begin());
	// epsilons.push_back(inf);
	while (vertices.size() < distances.size()) {
		for (Vertex v = distances.begin(); v != distances.end(); ++v)
			dist[v] = std::min(dist[v], distances(v, vertices.back()));
		auto max = std::max_element(dist.begin(), dist.end());
		vertices.push_back(max - dist.begin());
		epsilons.push_back(*max);
	}
	epsilons.push_back(0);
}

inline void compute_edges(const PairDistances& distances, 
				   const VertexVector& vertices, 
				   const EpsilonVector& epsilons,
				   const DistanceType& multiplier, 
				   EdgeVector& edges)
{
	for (unsigned i = 0; i != vertices.size(); ++i)
		for (unsigned j = i + 1; j != vertices.size(); ++j) {
			Vertex u = vertices[i];
			Vertex v = vertices[j];
			if (distances(u, v) <= multiplier * epsilons[j - 1]) edges.emplace_back(u, v);
		}
	std::sort(edges.begin(), edges.end(), 
		[&distances](const Edge& e1, const Edge& e2) {
			return distances(std::get<0>(e1), std::get<1>(e1)) < 
				   distances(std::get<0>(e2), std::get<1>(e2));
	});
}

inline void compute_positive_cofaces(
	const PairDistances& distances,
	const VertexVector& vertices, 
	const EdgeVector& edges, 
	const EpsilonVector& epsilons,
	const DistanceType& multiplier,
	Generator& rips,
	short unsigned& skeleton,
	unsigned& ce,
	unsigned& i,
	SimplexSet& cofaces)
{
	cofaces.clear();

	// Add anything else that needs to be inserted into the complex
	while (ce < edges.size()) {
		Vertex u, v;
		std::tie(u, v) = edges[ce];
		if (distances(u, v) <= multiplier * epsilons[i - 1])
			++ce;
		else
			break;
		// std::cout << "Adding cofaces of " << u << ' ' << v << std::endl;
		rips.edge_cofaces(
			u, v, 
			skeleton, 
			multiplier * epsilons[i - 1], 
			[&cofaces](Simplex&& s) { cofaces.insert(s); },
			vertices.begin(), 
			vertices.begin() + i + 1);
	}
}

inline void compute_negative_cofaces(
	const VertexVector& vertices, 
	const EpsilonVector& epsilons,
	const DistanceType& multiplier,
	Generator& rips,
	short unsigned& skeleton,
	unsigned& i,
	SimplexSet& cofaces)
{
	cofaces.clear();
	rips.vertex_cofaces(
		vertices[i], 
		skeleton, 
		multiplier * epsilons[i - 1], 
		[&cofaces](Simplex&& s) { cofaces.insert(s); },
		vertices.begin(), 
		vertices.begin() + i + 1);
	// std::cout << "Total cofaces: " << cofaces.size() << std::endl;
}

inline void build_rips_zigzag_filtration(std::vector<std::vector<int> > &simpls, 
								  std::vector<bool>& dirs, 
								  unsigned int numberOfPoints, 
								  int seed = -1,
								  short unsigned skeleton = 2, 
								  DistanceType multiplier = 6)
{
	// std::cout << "Building filtration" << std::endl;

	PointContainer 			points = compute_points(numberOfPoints, seed);

	// Construct distances and Rips generator
	PairDistances 			distances(points);
	Generator 				rips(distances);

	// Order vertices and epsilons (in maxmin fashion)
	VertexVector 	vertices;
	EpsilonVector 	epsilons;
	EdgeVector 		edges;
	
	compute_vertices_and_epsilons(distances, vertices, epsilons);

	// Generate and sort all the edges
	compute_edges(distances, vertices, epsilons, multiplier, edges);

	// Insert vertices
	for (auto v : vertices) {
		// Add a vertex
		simpls.push_back({static_cast<int>(v)});
		dirs.push_back(true);
	}

	// Process vertices
	dlog::progress	progress(vertices.size());
	unsigned 		ce = 0;		// index of the current one past last edge in the complex
	SimplexSet 		cofaces;	// record the cofaces of all the simplices that need to be removed and reinserted
	
	for (unsigned stage = 0; stage != vertices.size() - 1; ++stage) {
		unsigned i = vertices.size() - 1 - stage;

		/* Increase epsilon */
		compute_positive_cofaces(distances, vertices, edges, epsilons, multiplier, rips, skeleton, ce, i, cofaces);

		for (auto& s : cofaces) {
			// std::cout << "Inserting: " << s << std::endl;
			simpls.emplace_back(s.begin(), s.end());
			dirs.push_back(true);
		}

		/* Remove the vertex */
		// std::cout << "Removing vertex: " << vertices[i] << std::endl;
		compute_negative_cofaces(vertices, epsilons, multiplier, rips, skeleton, i, cofaces);

		for (auto& s : cofaces | ba::reversed) {
			// std::cout << "Removing: " << s << std::endl;
			simpls.emplace_back(s.begin(), s.end());
			dirs.push_back(false);
		}

		++progress;
	}

	// std::cout << "Finished" << std::endl;
}

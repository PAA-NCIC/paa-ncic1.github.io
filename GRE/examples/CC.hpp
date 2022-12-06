#ifndef _GRE_APP_CC_HPP
#define _GRE_APP_CC_HPP

#include "Types.hpp"
#include "VertexCode/VertexCode.hpp"
#include "Graph/VertexSet.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"

typedef vid_t vertex_data_t;
typedef vid_t scatter_data_t;
typedef vid_t combine_data_t;

typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef	GRE::VertexCode<vertex_data_t, scatter_data_t, combine_data_t, DistributedGraph_t> ConnectedComponents_Base;

//
//This algorithm is ONLY for undirected graphs.
//
class ConnectedComponents: public ConnectedComponents_Base{
public:
	ConnectedComponents(DistributedGraph_t& dg): ConnectedComponents_Base(dg), cc_scatter(NULL), cc(NULL){}
	~ConnectedComponents(){}
	void alloc(){
		allocScatterData(cc_scatter);
		allocCombineData(cc, std::numeric_limits<vid_t>::max());
		allocVertexData(cc);
	}
	void free(){
		freeScatterData(cc_scatter);
		freeCombineData(cc);
		freeVertexData(cc);
	}

public://Instantiation of Scatter-Combine Model
	//Init
	inline bool init(const lvid_t lv){
		cc[lv] = graph.get_vid(lv);
		cc_scatter[lv] = cc[lv];
		if(graph.getOutdegree(lv)>0){
			return true;
		} else {
			return false;
		}
	}
	
	//Scatter
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, cc_scatter[lsrc]);
	}

	//Combine
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		if(cc[ldst]>msg_data){
			lock.spinLock();
			if(cc[ldst]>msg_data){
				cc[ldst]=msg_data;
				lock.unLock();
				return true;
			}
			lock.unLock();
		}
		return false;
	}

	//Apply
	inline bool apply(const lvid_t lv){
		cc_scatter[lv] = cc[lv];
		return true;
	}
private:
	vertex_data_t* cc;
	scatter_data_t* cc_scatter;
};

#endif

#include <iostream>
#include <vector>
#include <string>

#include "Types.hpp"
#include "Graph/DistributedGraph.hpp"
#include "Graph/LocalGraph.hpp"
#include "Graph/GraphIngress.hpp"
#include "Graph/VertexSet.hpp"

#include "Comm/CommUnderlying.hpp"
#include "Comm/BufferingComm.hpp"

#include "Util/BufferPool.hpp"
#include "Engine/Engine.hpp"

//#include "App/BC.hpp"
//////////////////////////////////////////////


//vertex_data_t for BFS
struct vtxData{
	int level;
	int numPaths;	
	vtxData():level(-1), numPaths(0){}
	~vtxData(){}
};
typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef	GRE::VertexCode<int, int, vtxData, DistributedGraph_t> Forward_Base;
typedef	GRE::VertexCode<double, double, double, DistributedGraph_t> Backward_Base;

class Forward: public Forward_Base{
public:
	Forward(DistributedGraph_t& dg, vtxData*& _state, int _currLevel=1):Forward_Base(dg), state(_state), numPaths(NULL), sum(NULL),
		currLevel(_currLevel){}

	~Forward(){}

	void alloc(){
		allocScatterData(numPaths, 0);
		allocCombineData(sum, 0);
		allocVertexData(state);
	}
	void free(){
		freeScatterData(numPaths);
		freeCombineData(sum);
		//freeVertexData(state);
	}
	void aggregate(){
		currLevel++;
	}
	void setCurrLevel(int l){currLevel=l;}
	int getCurrLevel(){return currLevel;}
private:
	//vertex_data
	vtxData*& state;
	//scatter data
	int* numPaths;
	//combine data
	int* sum;
	//other
	int currLevel;
public:
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		//update my state
		state[lv].numPaths=1;
		state[lv].level=0;
		//update scatter data
		numPaths[lv]=state[lv].numPaths;
		return true;
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, numPaths[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		lock.spinLock();
		sum[ldst]+=msg_data;
		lock.unLock();
		return true;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		if(state[lv].level<0){//unvisited
			//update my state
			state[lv].level=currLevel;
			state[lv].numPaths=sum[lv];
			//update scatter data
			numPaths[lv]=state[lv].numPaths;
			return true;
		}
		return false;
	}

	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){
	}
	inline void reset_scatter(const lvid_t lsrc){
	}
}; 

class Backward: public Backward_Base{
public:
	Backward(DistributedGraph_t& dg, vtxData* _bfsState, double*& _dependencies, int _currLevel=0):
		Backward_Base(dg), bfsState(_bfsState), dependencies(_dependencies), 
		pathwgt(NULL), sum(NULL), currLevel(_currLevel){}
	~Backward(){}

	void setCurrLevel(const int l){ currLevel=l;}
	int getCurrLevel(){ return currLevel;}
	void aggregate(){}

	void alloc(){
		assert(bfsState!=NULL);
		allocScatterData(pathwgt, 0);
		allocCombineData(sum, 0);
		allocVertexData(dependencies, 0.0);
	}
	void free(){
		freeScatterData(pathwgt);
		freeCombineData(sum);
		//freeVertexData(dependencies);
	}
private:
	//external
	vtxData* bfsState;
	//vertex_data
	double*& dependencies;
	//scatter_data
	double* pathwgt;
	//combine_data
	double* sum;
	//other
	int currLevel;
public:
	///////////////////////////////Init//////////////////////////////////
	inline bool init(const lvid_t lv){
		if(bfsState[lv].level==currLevel){
			pathwgt[lv]=(1.0+dependencies[lv])/bfsState[lv].numPaths;
			return true;
		} else
			return false;
	}

	///////////////////////////////Scatter//////////////////////////////////
	template <typename Ctx_t, typename EdgeIterator_t>
	inline void scatter(Ctx_t& ctx, const lvid_t lsrc, const lvid_t ldst, const EdgeIterator_t& it){
		ctx.engine->sendMessage(ctx, ldst, pathwgt[lsrc]);
	}

	///////////////////////////////Combine//////////////////////////////////
	template <typename Lock_t>
	inline bool combine(const lvid_t ldst, const combine_data_t& msg_data, Lock_t& lock){
		lock.spinLock();
		sum[ldst]+=msg_data;
		lock.unLock();
		return true;
	}

	///////////////////////////////Apply//////////////////////////////////
	inline bool apply(const lvid_t lv){
		if(bfsState[lv].level==currLevel-1){
			dependencies[lv]=sum[lv]*bfsState[lv].numPaths;
			//pathwgt[lv]=(1.0+dependencies[lv])/numPaths[lv];
			//return true;
		}
		return false;
	}
	///////////////////////////////Misc//////////////////////////////////
	inline void reset_combiner(const lvid_t ldst){
		sum[ldst]=0.0;
	}
	inline void reset_scatter(const lvid_t lsrc){
	}
}; 


//////////////////////////////////////////////
typedef GRE::Util::BufferPool BufferPool_t;
//typedef GRE::Graph::DistributedGraph<GRE::Graph::GraphInCSR> DistributedGraph_t;
typedef GRE::Graph::SourceVertexSet SourceVertexSet_t;

typedef Forward Forward_t;
typedef Backward Backward_t;
typedef GRE::SynchronousEngine<Forward_t> ForwardEngine_t;
typedef GRE::SynchronousEngine<Backward_t> BackwardEngine_t;

int main(int argc, char** argv)
{	
	GRE::COMM::Underlying comm;
	comm.init(argc, argv);

	BufferPool_t pool;
	pool.config(2048, 8192*256, true);

	std::string prefix(argv[1]);

	///////////////////////////
	//Load graph
	DistributedGraph_t distributedGraph;
	comm.barrier();
	distributedGraph.ingress<InOut>(prefix, comm, pool, "ascii", "edgelists");
	comm.barrier();
	//load source vertices
	SourceVertexSet_t srcSet;
	std::vector<vid_t> src;
	const int n = GRE::Util::loadFromFile(prefix+".roots", src, 0);
	if(n==0){
		std::cerr<<"Err in reading src vertex from file."<<std::endl;
		exit(0);
	}
	srcSet.addVertex(src[0], distributedGraph);
	///////////////////////////

	///////////////////////////
	//phase-1:forward bfs
	vtxData* vtxState=NULL;//allocated in bfs
	Forward_t bfs(distributedGraph, vtxState);
	ForwardEngine_t bfsengine(distributedGraph, comm, pool, bfs);
	bfsengine.initialize();
	bfsengine.run<Out>(srcSet, GRE::Par_ThreadPool);
	bfsengine.output_stat();
	comm.outputProfilingData();
	bfsengine.finalize();
	///////////////////////////

	///////////////////////////
	int nLevels=bfs.getCurrLevel();
	///////////////////////////

	///////////////////////////
	//phase-2:backward traces
	double* dependencies=NULL;
	Backward_t bp(distributedGraph, vtxState, dependencies);
	BackwardEngine_t bpengine(distributedGraph, comm, pool, bp);
	bpengine.initialize();
	for(int i=nLevels-1; i>=0; i--){
		srcSet.setFULL();
		bp.setCurrLevel(i);
		bpengine.run<In>(srcSet, GRE::Par_ThreadGroup);
	}
	bpengine.finalize();
	///////////////////////////

	///////////////////////////
	//dependencies[] should be save to file.
	delete []vtxState;
	delete []dependencies;

	///////////////////////////
	comm.finalize();
	pool.finalize();
}

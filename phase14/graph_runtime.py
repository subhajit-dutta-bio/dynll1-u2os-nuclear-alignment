"""Optional launch-overhead optimization of the unchanged network operations."""
import torch
def cache_fixed_relative_positions():
    # These values depend only on frozen model parameters and tile dimensions.
    # Keeping them on the GPU avoids CPU index transfers inside graph capture.
    import segment_anything.modeling.image_encoder as encoder
    previous=encoder.get_rel_pos;cache={}
    def cached(q_size,k_size,rel_pos):
        key=(q_size,k_size,rel_pos.data_ptr(),rel_pos.dtype,rel_pos.device,rel_pos._version)
        if key not in cache:cache[key]=previous(q_size,k_size,rel_pos)
        return cache[key]
    encoder.get_rel_pos=cached
    return cache
class GraphForward:
    def __init__(self,net):
        self.native=net.forward;self.cache={};self.pool=torch.cuda.graph_pool_handle()
    def __call__(self,x):
        key=(tuple(x.shape),x.dtype,x.device)
        if key not in self.cache:
            static=torch.zeros_like(x);stream=torch.cuda.Stream();stream.wait_stream(torch.cuda.current_stream())
            with torch.cuda.stream(stream),torch.no_grad():
                for _ in range(3):self.native(static)
            torch.cuda.current_stream().wait_stream(stream)
            graph=torch.cuda.CUDAGraph()
            with torch.cuda.graph(graph,pool=self.pool),torch.no_grad():out=self.native(static)
            self.cache[key]=(graph,static,out)
        graph,static,out=self.cache[key];static.copy_(x);graph.replay();return out
def install(net):
    positions=cache_fixed_relative_positions();wrapper=GraphForward(net);wrapper.relative_positions=positions;net.forward=wrapper;return wrapper

from pipelines.base_layer import Layer, pipe_object

class Input(object):
    # Prepare inputs pipe object for a layer using normal keywords
    def __init__(self, layer_name="get_count_from_sra_list", **kwargs ):
        self.layer_name = layer_name

    def __call__(self, **kwargs):
        out = vars(self.__class__)[self.layer_name](self,**kwargs)
        return out

    def check_kwargs_error(self, kwargs, allowed_kwargs):
        for kwarg in kwargs.keys():
            if kwarg not in allowed_kwargs:
                raise TypeError('Keyword argument not understood:{}. \nAllowed kwargs:{}'.
                                format(kwarg, allowed_kwargs))
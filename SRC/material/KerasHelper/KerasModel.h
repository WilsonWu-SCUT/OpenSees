#pragma once

namespace keras2cpp
{
    class KerasModel 
    {
    public:
        KerasModel();
        KerasModel(const std::string& file_path);
        ~KerasModel();

    public:
        std::vector<float> predict(const std::vector<float>& input_param_vec);

    private:
        //keras模型指针
        Model* keras_model_ptr_;
    };
}
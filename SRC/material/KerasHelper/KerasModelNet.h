#pragma once

namespace Keras
{
    public ref class KerasModelNet 
    {
    public:
        KerasModelNet();
        KerasModelNet(String^ file_path);
        ~KerasModelNet();

    public:
        List<double>^ Predict(List<double>^ inputParamList);

    private:
        //keras模型指针
        KerasModel* keras_model_ptr_;
    };
}
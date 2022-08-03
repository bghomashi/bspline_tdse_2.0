#ifndef __FACTORY_H__
#define __FACTORY_H__

// *************************************************************
// In addition to virtual functions (Startup, Shutdown, Compute)
//    the final derived class to be created by the Factory must
//    implement three(3) static functions with the signatures:
//
//      static std::shared_ptr<T> Create(const nlohmann::json& )
//      static bool Validate(const nlohmann::json& )
//      static std::string GetName();
//
//    Create - the function the factory calls to create the object
//    from an input json file.
//    Validate - validates the input json file.
//    GetName - return the name the factory will use to reference
//    this object
// ************************************************************

#include "common/utility/json.hpp"
#include <unordered_map>
#include <string>
#include <functional>
#include <memory>

namespace utl {
    template <typename T>
    struct Factory
    {
        typedef std::shared_ptr<T> ObjectPtr_t;
        typedef std::function<ObjectPtr_t(const nlohmann::json &)> CreateFunc_t;
        typedef std::function<bool(const nlohmann::json &)> ValidateFunc_t;
        typedef std::unordered_map<std::string, CreateFunc_t> CreateMap_t;
        typedef std::unordered_map<std::string, ValidateFunc_t> ValidateMap_t;

        template <typename S>
        struct Register : public T
        {
            Register() { (void)s_register; }

            static bool registerT()
            {
                Factory::AddCreateFunc(S::GetName(), &S::Create);
                Factory::AddValidateFunc(S::GetName(), &S::Validate);
                return true;
            }

        protected:
            static bool s_register;
        };

        static ObjectPtr_t Create(const std::string &name, const nlohmann::json &input);
        static bool Validate(const std::string &name, const nlohmann::json &input);
        static bool Exists(const std::string &name);
        static bool AddCreateFunc(const std::string &name, CreateFunc_t f);
        static bool AddValidateFunc(const std::string &name, ValidateFunc_t f);

    private:
        static ValidateMap_t &ValidateFunctions();
        static CreateMap_t &CreateFunctions();
    };

    template <typename T>
    template <typename S>
    bool Factory<T>::Register<S>::s_register = Factory<T>::Register<S>::registerT();

    // ----- function definitions --------
    template <typename T>
    typename Factory<T>::ObjectPtr_t Factory<T>::Create(const std::string &name, const nlohmann::json &input)
    {
        auto it = CreateFunctions().find(name);
        if (it != CreateFunctions().end())
            return it->second(input);
        return nullptr;
    }
    template <typename T>
    bool Factory<T>::Validate(const std::string &name, const nlohmann::json &input)
    {
        auto it = ValidateFunctions().find(name);
        if (it != ValidateFunctions().end())
            return it->second(input);
        return false;
    }
    template <typename T>
    bool Factory<T>::Exists(const std::string &name)
    {
        return (CreateFunctions().find(name) != CreateFunctions().end());
    }
    template <typename T>
    bool Factory<T>::AddCreateFunc(const std::string &name, CreateFunc_t f)
    {
        CreateFunctions().insert(std::make_pair(name, f));
        return true;
    }
    template <typename T>
    bool Factory<T>::AddValidateFunc(const std::string &name, ValidateFunc_t f)
    {
        ValidateFunctions().insert(std::make_pair(name, f));
        return true;
    }

    template <typename T>
    typename Factory<T>::ValidateMap_t &Factory<T>::ValidateFunctions()
    {
        static ValidateMap_t s_validateFunctions;
        return s_validateFunctions;
    }
    template <typename T>
    typename Factory<T>::CreateMap_t &Factory<T>::CreateFunctions()
    {
        static CreateMap_t s_createFunctions;
        return s_createFunctions;
    }
}

#endif

#include "appmanager.h"

#include <iostream>
#include <chrono>
#include <thread>

AppManager::AppManager()
{
    std::cout << __PRETTY_FUNCTION__ << " created " << std::endl;
}

AppManager::~AppManager()
{
    clean_up();
}

int AppManager::run() noexcept(false)
{
    std::cout << "Operating..." << std::endl;
    
    std::this_thread::sleep_for(std::chrono::seconds(2));

    if (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) < 0.1f)
        throw std::runtime_error("Example error during runtime");

    std::cout << "All good!" << std::endl;


    if (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) < 0.1f)
        return 1;

    return 0;
}

void AppManager::clean_up()
{
    std::cout << __PRETTY_FUNCTION__ << " is being destroyed " << std::endl;
}
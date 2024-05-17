#ifndef APPMANAGER_H
#define APPMANAGER_H

#include <stdexcept>

class AppManager
{
public:
    AppManager();
    ~AppManager();

    int run() noexcept(false);
private:

};

#endif // APPMANAGER_H
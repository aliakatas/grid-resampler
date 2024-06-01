#ifndef APPMANAGER_H
#define APPMANAGER_H

#include <stdexcept>

class AppManager
{
public:
    AppManager();
    ~AppManager();

    int run() noexcept(false);

    void clean_up();
private:

};

#endif // APPMANAGER_H
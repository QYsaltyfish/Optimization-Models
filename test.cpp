//
// Created by QY on 2024/10/23.
//

// You may need to build the project (run Qt uic code generator) to get "ui_test.h" resolved

#include "test.h"
#include "ui_test.h"


test::test(QWidget *parent) :
        QWidget(parent), ui(new Ui::test) {
    ui->setupUi(this);
}

test::~test() {
    delete ui;
}
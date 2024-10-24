//
// Created by QY.
//

#include <QApplication>
#include <QPushButton>
#include <QVBoxLayout>
#include <QWidget>
#include <QComboBox>
#include <QFileDialog>
#include <QLabel>
#include "SimplexAlgorithm.h"

class SolverApp : public QWidget {
public:
    SolverApp() {
        // Create layout
        auto *layout = new QVBoxLayout(this);

        // Create problem type selection box
        problemTypeCombo = new QComboBox();
        problemTypeCombo->addItems({"Linear Programming", "Integer Programming"});

        // Create solution method selection box
        methodCombo = new QComboBox();
        updateMethodCombo();  // Initialize methodCombo based on the default selection

        // Connect the signal to the slot
        connect(problemTypeCombo, &QComboBox::currentIndexChanged, this,
                &SolverApp::updateMethodCombo);

        layout->addWidget(new QLabel("Select Problem Type:"));
        layout->addWidget(problemTypeCombo);
        layout->addWidget(new QLabel("Select Solution Method:"));
        layout->addWidget(methodCombo);

        // Create button to select file
        fileButton = new QPushButton("Select CSV File");
        layout->addWidget(fileButton);
        connect(fileButton, &QPushButton::clicked, this, &SolverApp::selectFile);

        // Create solve button
        solveButton = new QPushButton("Solve");
        layout->addWidget(solveButton);
        connect(solveButton, &QPushButton::clicked, this, &SolverApp::solve);

        // Set window properties
        setLayout(layout);
        setWindowTitle("Solver");
        resize(300, 200);
    }

private slots:
    void selectFile() {
        QString fileName = QFileDialog::getOpenFileName(this, "Select CSV File",
                                                        "", "CSV Files (*.csv)");
        if (!fileName.isEmpty()) {
            // Store the file path
            csvFilePath = fileName.toStdString();
        }
    }

    void solve() {
        QString problemType = problemTypeCombo->currentText();
        QString method = methodCombo->currentText();

        // Create a label to display the output
        if (OutputLabel == nullptr) {
            OutputLabel = new QLabel();
            layout()->addWidget(OutputLabel);
        }

        // Call the solve function with problem type, method, and file path
        if (problemType == "Linear Programming") {
            if (method == "Simplex Method") {
                simplexMethod();
            }
        }
    }

private:
    QComboBox *problemTypeCombo;
    QComboBox *methodCombo;
    QPushButton *fileButton;
    QPushButton *solveButton;
    QLabel *OutputLabel = nullptr;
    std::string csvFilePath;  // Store the selected CSV file path

    void updateMethodCombo() {
        methodCombo->clear();  // Clear current items
        QString problemType = problemTypeCombo->currentText();

        if (problemType == "Linear Programming") {
            methodCombo->addItems({"Simplex Method", "Interior Point Method (not yet finished)"});
        } else if (problemType == "Integer Programming") {
            methodCombo->addItems({"Not yet finished."});
        }
    }

    void simplexMethod() {
        try {
            LpProblem problem = LpProblem(csvFilePath);
            OutputLabel->setText("Solving...");
            auto simplex = SimplexAlgorithm(problem);
            double output = simplex.solve();
            if (output != inf) {
                OutputLabel->setText(QString("Optimal Value: %1").arg(output));
            } else {
                OutputLabel->setText("The problem is unbounded");
            }
        } catch (const std::invalid_argument& e) {
            OutputLabel->setText("Error: " + QString::fromStdString(e.what()));
        }
    }
};

int main(int argc, char *argv[]) {
    QApplication a(argc, argv);
    SolverApp app;
    app.show();
    return QApplication::exec();
}
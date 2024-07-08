# Load necessary libraries
library(tidyverse)    # For data manipulation and visualization
library(caret)        # For model training and evaluation
library(glmnet)       # For elastic net regularization
library(broom)        # For tidying up model outputs
library(pROC)         # For ROC and AUC

# Step 1: Load your data
# Replace 'your_data.csv' with the path to your dataset
data <- read.csv('your_data.csv')

# Step 2: Explore the data
# Display the first few rows of the dataset
head(data)

# Summary statistics of the data
summary(data)

# Check for missing values
colSums(is.na(data))

# Step 3: Data Preprocessing
# Convert categorical variables to factors if needed
# Example: data$Category <- as.factor(data$Category)

# Remove or impute missing values if any
# Example: data <- na.omit(data)

# Define predictor and response variables
# Assuming the target variable is named 'TargetVariable' and predictors are all other columns
response <- data$TargetVariable
predictors <- data %>% select(-TargetVariable)

# Split the data into training and testing sets
set.seed(123)  # For reproducibility
trainIndex <- createDataPartition(response, p = 0.8, list = FALSE)
trainData <- predictors[trainIndex, ]
trainResponse <- response[trainIndex]
testData <- predictors[-trainIndex, ]
testResponse <- response[-trainIndex]

# Step 4: Fit the elastic net model with cross-validation
# Convert data to matrix format required by glmnet
x_train <- as.matrix(trainData)
y_train <- as.factor(trainResponse)
x_test <- as.matrix(testData)

# Fit the model using cross-validation to find the best lambda
cv_model <- cv.glmnet(x_train, y_train, alpha = 0.5, family = "multinomial")

# Step 5: Model summary
# Best lambda value
best_lambda <- cv_model$lambda.min
print(paste("Best lambda:", best_lambda))

# Coefficients of the model at the best lambda
coef(cv_model, s = "lambda.min")

# Step 6: Make predictions on the test set
predictions <- predict(cv_model, newx = x_test, s = "lambda.min", type = "response")

# Convert probabilities to class labels
predicted_classes <- apply(predictions, 1, which.max)
predicted_classes <- levels(y_train)[predicted_classes]

# Step 7: Evaluate the model
# Confusion matrix
confusionMatrix(factor(predicted_classes), factor(testResponse))

# ROC curve and AUC (for each class)
for (i in 1:length(levels(y_train))) {
  roc_curve <- roc(as.numeric(testResponse == levels(y_train)[i]), predictions[, i])
  plot(roc_curve, main = paste("ROC Curve for class", levels(y_train)[i]))
  print(paste("AUC for class", levels(y_train)[i], ":", auc(roc_curve)))
}

# Step 8: Tidy model output (optional)
tidy_model <- tidy(cv_model)
glance(cv_model)
augment(cv_model, newdata = x_test)

# Save the model if needed
# saveRDS(cv_model, file = "elastic_net_model.rds")

# Load the model later if needed
# loaded_model <- readRDS("elastic_net_model.rds")

%%  INITIATING INPUTS! 

inputs.nmin1 = 3; 
inputs.nmin2 = 1; 
inputs.OnlyTree = 1; 
inputs.Tria = 0; 
inputs.Dist = 1; 
inputs.MinCylRad = 0.0025;
inputs.ParentCor = 1; 
inputs.TaperCor = 1; 
inputs.GrowthVolCor = 0; 
inputs.GrowthVolFac = 1.5; 
inputs.filter.k = 10;
inputs.filter.radius = 0.00;
inputs.filter.nsigma = 1.5;
inputs.filter.PatchDiam1 = 0.05;
inputs.filter.BallRad1 = 0.075;
inputs.filter.ncomp = 2;
inputs.filter.EdgeLength = 0.004;
inputs.filter.plot = 1;
inputs.name = 'tree'; 
inputs.tree = 1;
inputs.model = 1;
inputs.savemat = 1;
inputs.savetxt = 1; 
inputs.plot = 0;
inputs.disp = 0; 



% User Defined 
% -------------

train_set = 34; % ~ 70%
test_set  = 15; % ~ 30%
no_of_iterations = 5;   % Number TRIALS OR ITERATIONS!
no_of_model_runs = 5;   % Each RUN for each tree will be repeated 5 times and the average DBH will be taken, The reason is that TreemQSM gives different results even for the same inputs.

DBH_filepath = 'OriID_vs_TreeQSMID.xlsx'; % YOU NEED TO CHANGE THIS! Directory where DBH observed Data (Excel sheet) is located.

train_range = ['E' num2str(2) ':E' num2str(train_set+1)]; % I am specifying the cells where the DBH data is stored ((from D2 to D21))
test_range  = ['E' num2str(train_set+2) ':E' num2str(train_set+test_set+1)]; % Extract the cells containing test data.

disp("------------------------")
disp("Define Inputs Completed!")
disp("------------------------")

%%
lasDirectory = ''; % This means that the .las files are in the main src directory.
lasFiles = dir(fullfile(lasDirectory, '*.las'));
file_nums = cellfun(@(x) str2double(x(1:end-4)), {lasFiles.name});
[~, sorted_indices] = sort(file_nums);
lasFiles = lasFiles(sorted_indices);
variable_names = cell(1, train_set); % This is just to assign variable names to each point cloud before saving. 

for tree = 1:train_set
    name = lasFiles(tree).name;
    lasReader = lasFileReader(name);
    ptCloud = readPointCloud(lasReader);
    P = ptCloud.Location - mean(ptCloud.Location);
    var_name = ['P_' num2str(tree)];
    eval([var_name ' = P;']);
    variable_names{tree} = var_name;
end

save('trees_train', variable_names{:});

disp("-----------------------")
disp("Saving Trees Completed!")
disp("-----------------------")

%%

rmse_old=1e6; % Initiating a large value for RMSE just for the sake of comparison and storing the smallest RMSE. 
rmse_for_iteration= zeros(no_of_iterations, 1); % Creating an empty array to store RMSE for each trial (for monitoring purposes).

DBH_obs = xlsread(DBH_filepath, train_range); % Here We're importing the array of the observed DBH.From Sheet 2. 

for iteration = 1:no_of_iterations % START THE TRIALS LOOP.

    inputs.PatchDiam1    = 0.05 + (0.10 - 0.05) * rand;  % The 0.05 is the lower limit and the 0.10 is the upper limit. [0.05 - 0.10]
    inputs.PatchDiam2Min = 0.02 + (0.06 - 0.02) * rand;  % The 0.02 is the lower limit and the 0.06 is the upper limit. [0.02 - 0.06]
    inputs.PatchDiam2Max = 0.03 + (0.15 - 0.03) * rand;  % The 0.03 is the lower limit and the 0.15 is the upper limit. [0.03 - 0.15]
    inputs.BallRad1 = inputs.PatchDiam1 + 0.015;
    inputs.BallRad2 = inputs.PatchDiam2Max + 0.01;

    QSMs = make_models_parallel( 'trees_train' , 'QSMs trees' , no_of_model_runs , inputs );

    % Since the output data structure is complicated, I have to use nested loop to get the output DBH, the below nested loop is just averaging the DBH.
    DBH_qsm=zeros(train_set,1);
    final_params=zeros(4,1); % Because we are dealing with 3 parameters that we need to store and monitor and the fourth is the RMSE score 

    kk=1;
    for i = 1 : no_of_model_runs : no_of_model_runs*train_set
        DBH_array=[];

        c=1;
            
            for j = i:i+no_of_model_runs-1
                DBH = double(QSMs(j).treedata.DBHqsm);
                DBH_array(c)= DBH;
            end
            
            avg_DBH=mean(DBH_array);

        DBH_qsm(kk)=avg_DBH;
        kk = kk+1;
    end

    rmse_score = sqrt(mean((DBH_obs - DBH_qsm).^2));
    disp(rmse_score)
    rmse_for_iteration(iteration)=rmse_score; % HERE I am storing the RMSE for each trial over the 49 trees (training set).
    if rmse_score < rmse_old
        
        % Store values in a 4x1 vector
        final_params = [inputs.PatchDiam1; inputs.PatchDiam2Min; inputs.PatchDiam2Max; rmse_score]; % FINAL RESULTS CONVENTION.
        rmse_old = rmse_score;

    end % End for the if statement.
    disp(['--------Iteration (' num2str(iteration) ') Completed']);
end


%%

variable_names_t = cell(1, test_set); % This is just to assign variable names to each point cloud before saving. 

cc=1;
for tree = train_set+1:train_set+test_set
    name_t = lasFiles(tree).name;
    lasReader_t = lasFileReader(name_t);
    ptCloud_t = readPointCloud(lasReader_t);
    P_t = ptCloud_t.Location - mean(ptCloud_t.Location);
    var_name_t = ['Pt_' num2str(tree)];
    eval([var_name_t ' = P_t;']);
    variable_names_t{cc} = var_name_t;
    cc = cc+1;
end

save('trees_test', variable_names_t{:});


inputs.PatchDiam1    = final_params(1); % NOW 
inputs.PatchDiam2Min = final_params(2); % THEY
inputs.PatchDiam2Max = final_params(3); % ARE FIXED! 
inputs.BallRad1 = inputs.PatchDiam1+0.015;
inputs.BallRad2 = inputs.PatchDiam2Max+0.01;

disp(["Test DBH From Cells:" test_range])
DBH_obs_test = xlsread(DBH_filepath, test_range); % Here We're importing the array of the observed DBH But now for testing.

QSMs_test = make_models_parallel( 'trees_test' , 'QSMs trees test' , no_of_model_runs , inputs );

k=1;
DBH_qsm_test=[];
for i = 1 : no_of_model_runs : no_of_model_runs*test_set
DBH_array=[];
    c=1;
    for j = i:i+no_of_model_runs-1
        DBH = double(QSMs(j).treedata.DBHqsm);
        DBH_array(c)= DBH;
    c=c+1;
    end
    avg_DBH=mean(DBH_array);

DBH_qsm_test(k)=avg_DBH;
k = k+1;
end

final_test_result(:,1)= 100*DBH_qsm_test;
final_test_result(:,2)= DBH_obs_test;
disp("Final Testing Results:")
disp("   -----------------")
disp("   DBH_qsm   DBH_obs")
disp("   -------   -------")
disp(final_test_result)

%% Midul HW 9 %%
%%%%%%%%%%%%%%%%
% This is the main script for HW 9.
% To see the written answers and explanations, see PDF in the same zip file.

clear all; close all; clc;

%% Q1 - Hopfield

% Sections a - d
N = 120; P = [3, 7, 25]; lambda = 0; tau = 0;
iterations = 100; repeats = 100; 
noise = 0; temporal_association = 0;

figure;
for i = 1:length(P)
    subplot(length(P),1,i)
    stability_check = stability(N, P(i), iterations, repeats, noise, temporal_association, tau, lambda);
    plot(1:iterations, stability_check); xlabel('Repeats'); ylabel('Overlap');
    title(['Overlap/Iterations for ',num2str(P(i)),' Initial Memories'])
    percent_1(i) = sum(stability_check == 1)/iterations;
end

% Change conditions - sections e
noise = 1; repeats = 1;

figure;
for i = 1:length(P)
    subplot(length(P),1,i)
    stability_check = stability(N, P(i), iterations, repeats, noise, temporal_association, tau, lambda);
    plot(1:iterations, stability_check); xlabel('Iterations'); ylabel('Overlap');
    title(['Overlap/Iterations for ',num2str(P(i)),' Initial memories, Noisy']);
    percent_2(i) = sum(stability_check == 1)/iterations;
end

%% Q2 -  Temporal Association

P = 10; iterations = 200;
lambda = 1.3; tau = 19; repeats = 1; 
noise = 0; temporal_association = 1;
connections_mat = stability(N, P, iterations, repeats, noise, temporal_association, tau, lambda);

memory_num = cell(1,P);
for k = 1:length(memory_num)
    memory_num{k} = ['Memory No.', num2str(k-1)];
end

figure;
hold on;
for p = 1:length(P)
    plot(1:iterations, connections_mat(:,:))
end
hold off; title(['Overlap/Iterations, lambda =',num2str(lambda),', tau =',num2str(tau)]);
xlabel('Iterations');ylabel('Overlap'); legend(memory_num, "Location", "northeastoutside");

%% Local Functions

% Checking system stability
function result = stability(N, P, iterations, repeats, noise, temporal_association, tau, lambda)

    if temporal_association == 1
        result = zeros(P,iterations);
    else
        if repeats == 1
            result = zeros(1, iterations);
        else
            result = zeros(1, repeats);
        end
    end

    for repeats = 1:repeats  
        S_mat = 2*round(rand(N,P)) - 1;
        state = S_mat(:,1); J = zeros(N);

        for j = 1:length(J)
           for n = 1:length(J)
              if n ~= j
                J(j,n) = (S_mat(j,:)*S_mat(n,:)')/N;
              end
           end
        end
        
        % Noise
        if noise == 1
            for j = 1:N
                rand_noise = randi(10);
                if rand_noise == 1
                    state(j) = (-1)*S_mat(j,1);
                end
            end
        end
        
        % Temporal Association
        if temporal_association == 1
            S_s = zeros(N,iterations);
            J_j = zeros(N);
            for j = 1:length(J_j)
                for n = 1:length(J_j) 
                    if n ~= j
                        J_j(j,n) = (S_mat(j,2:end)*S_mat(n,1:end-1)')/N;
                    end
                end
            end
        end

        for n = 1:iterations
            for k = 1:N
                if (temporal_association == 1) && (n > tau)
                    h = lambda*(J_j(k,:)*S_s(:,n - tau));
                else
                    h = 0;
                end
                h = h + J(k,:)*state; state(k) = sign(h);
            end

            if repeats == 1
                if temporal_association == 1
                    S_s(:,n) = state;
                    result(:,n) = (state(:,1)'*S_mat)'/N;
                else
                    result(n) = (state'*S_mat(:,1))/N;
                end
            end

        end

        if repeats ~= 1
            result(repeats) = (state'*S_mat(:,1))/N;
        end
    end
end